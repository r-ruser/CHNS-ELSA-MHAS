options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(WeightIt);library(survey);library(sandwich)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal");for(z in c(dd,od,td,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"03c_longitudinal_overlap_MSM.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
long<-readRDS(file.path(dd,"causal_treatment_long_internal_restricted.rds"))%>%arrange(Idind,wave)
oldw<-readRDS(file.path(dd,"causal_participant_weights_internal_restricted.rds"))

# Generalized overlap tilting for the observed three-state treatment. For state
# j with generalized propensity e_j(L), WeightIt ATO uses the multi-treatment
# overlap tilt h(L)={sum_k 1/e_k(L)}^{-1} and w_j=h/e_j. Products of normalized
# decision-specific weights define the longitudinal overlap population.
form<-fuel_state~previous_fuel+age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
form0<-fuel_state~age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
out<-list();audit<-list()
for(yr in sort(unique(long$wave))){
 z<-long%>%filter(wave==yr)%>%droplevels(); f<-if(all(z$previous_fuel=="entry"))form0 else form
 ow<-WeightIt::weightit(f,data=z,method="glm",estimand="ATO")
 z$overlap_wave_weight<-ow$weights/mean(ow$weights)
 if(any(!is.finite(z$overlap_wave_weight))||any(z$overlap_wave_weight<0))stop("Invalid overlap weight at ",yr)
 out[[as.character(yr)]]<-z
 audit[[as.character(yr)]]<-data.frame(wave=yr,N=nrow(z),states=n_distinct(z$fuel_state),mean=mean(z$overlap_wave_weight),SD=sd(z$overlap_wave_weight),minimum=min(z$overlap_wave_weight),P1=quantile(z$overlap_wave_weight,.01),P99=quantile(z$overlap_wave_weight,.99),maximum=max(z$overlap_wave_weight),ESS=sum(z$overlap_wave_weight)^2/sum(z$overlap_wave_weight^2))
}
owlong<-bind_rows(out)%>%arrange(Idind,wave)%>%group_by(Idind)%>%mutate(overlap_treatment_weight_cum=cumprod(overlap_wave_weight))%>%ungroup()
owfinal<-owlong%>%group_by(Idind)%>%summarise(overlap_treatment_weight=last(overlap_treatment_weight_cum),commid=last(commid),HD_z=last(HD_z),KDM_BAA=last(KDM_BAA),.groups="drop")%>%left_join(oldw%>%select(Idind,censoring_weight,strategy),by="Idind")%>%mutate(overlap_combined_weight=overlap_treatment_weight*censoring_weight)
saveRDS(owlong,file.path(dd,"causal_overlap_long_internal_restricted.rds"),compress="gzip");saveRDS(owfinal,file.path(dd,"causal_overlap_participant_weights_internal_restricted.rds"),compress="gzip")
write.csv(bind_rows(audit),file.path(od,"overlap_wave_weight_audit.csv"),row.names=FALSE)

wmean<-function(x,w)sum(x*w,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)
wvar<-function(x,w){m<-wmean(x,w);sum(w*(x-m)^2,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)}
pair<-function(x,g,w,a,b){i<-g==a&is.finite(x);j<-g==b&is.finite(x);den<-sqrt((wvar(x[i],w[i])+wvar(x[j],w[j]))/2);if(!is.finite(den)||den==0)return(NA_real_);(wmean(x[j],w[j])-wmean(x[i],w[i]))/den}
maxsmd<-function(x,g,w){lev<-unique(as.character(g));ps<-combn(lev,2,simplify=FALSE);v<-vapply(ps,function(p)pair(x,as.character(g),w,p[1],p[2]),numeric(1));if(all(!is.finite(v)))NA_real_ else max(abs(v[is.finite(v)]))}
vars<-c("age","urban_i","education_i","income_i","assets_i","smoking_i","alcohol_i")
bal<-list()
for(yr in sort(unique(owlong$wave))){z<-owlong%>%filter(wave==yr);for(v in vars)bal[[length(bal)+1]]<-data.frame(wave=yr,variable=v,SMD_unweighted=maxsmd(z[[v]],z$fuel_state,rep(1,nrow(z))),SMD_overlap=maxsmd(z[[v]],z$fuel_state,z$overlap_treatment_weight_cum))}
balance<-bind_rows(bal)%>%mutate(abs_overlap=abs(SMD_overlap),meets_0_10=abs_overlap<.10)
bsum<-balance%>%group_by(wave)%>%summarise(max_abs_SMD=max(abs_overlap,na.rm=TRUE),n_over_0_10=sum(abs_overlap>=.10,na.rm=TRUE),.groups="drop")
write.csv(balance,file.path(td,"Table_causal_balance_overlap.csv"),row.names=FALSE);write.csv(bsum,file.path(od,"overlap_balance_wave_summary.csv"),row.names=FALSE)

primary<-owfinal%>%filter(strategy%in%c("continued_any_solid","sustained_clean_only"))%>%mutate(A=as.integer(strategy=="sustained_clean_only"))
# Overlap treatment weights are intrinsically bounded per decision; IPCW is
# kept untruncated first, with 1/99 combined truncation as a prespecified check.
qq<-quantile(primary$overlap_combined_weight,c(.01,.99));primary<-primary%>%mutate(overlap_combined_weight_99=pmin(pmax(overlap_combined_weight,qq[1]),qq[2]))
ess<-primary%>%group_by(strategy)%>%summarise(N=n(),ESS_treatment=sum(overlap_treatment_weight)^2/sum(overlap_treatment_weight^2),ESS_combined99=sum(overlap_combined_weight_99)^2/sum(overlap_combined_weight_99^2),.groups="drop")
write.csv(ess,file.path(od,"overlap_primary_group_ESS.csv"),row.names=FALSE)
pass_balance<-all(bsum$max_abs_SMD<.10,na.rm=TRUE);pass_ess<-all(ess$ESS_combined99>=50);passed<-pass_balance&&pass_ess

estimand<-data.frame(name="Longitudinal generalized overlap average treatment effect",population="Eligible linked participants with high joint probability of receiving each observed fuel state at every treatment decision",contrast="Transition to and sustain clean-only versus continued any-solid",notation="E_OW[Y^B]-E_OW[Y^A]",difference_from_original="Changes the target population from the full eligible linked cohort to the longitudinal overlap population; it is not the original ATE.")
write.csv(estimand,file.path(td,"Table_overlap_estimand.csv"),row.names=FALSE)

if(passed){
 fitone<-function(y){des<-svydesign(ids=~commid,weights=~overlap_combined_weight_99,data=primary,nest=TRUE);fit<-svyglm(reformulate("A",response=y),design=des);b<-coef(fit)["A"];se<-sqrt(vcov(fit)["A","A"]);data.frame(outcome=y,estimate=b,lower95=b-qnorm(.975)*se,upper95=b+qnorm(.975)*se,p_value=2*pnorm(-abs(b/se)),N=nrow(primary),clusters=n_distinct(primary$commid))}
 msm<-bind_rows(fitone("HD_z"),fitone("KDM_BAA"));write.csv(msm,file.path(td,"Table3_overlap_MSM_effects.csv"),row.names=FALSE)
}else{
 msm<-data.frame(outcome=c("HD_z","KDM_BAA"),estimate=NA_real_,lower95=NA_real_,upper95=NA_real_,p_value=NA_real_,N=nrow(primary),clusters=n_distinct(primary$commid),reason="Overlap weighting failed prespecified balance or ESS threshold; outcome model withheld")
 write.csv(msm,file.path(td,"Table3_overlap_MSM_effects.csv"),row.names=FALSE)
}
status<-data.frame(method="longitudinal generalized overlap weighting",changed_estimand=TRUE,balance_pass=pass_balance,ESS_pass=pass_ess,overall_pass=passed,decision=ifelse(passed,"Overlap-population MSM fitted","Outcome model withheld; retain causal analysis as unstable"))
write.csv(status,file.path(od,"overlap_final_status.csv"),row.names=FALSE)
print(bsum);print(ess);print(status);if(passed)print(msm);cat("LONGITUDINAL_OVERLAP=",ifelse(passed,"PASS","FAIL"),"\n",sep="")
