options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(WeightIt);library(nnet);library(splines)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal");for(z in c(dd,od,td,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"03b_longitudinal_CBPS.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
long<-readRDS(file.path(dd,"causal_treatment_long_internal_restricted.rds"))%>%arrange(Idind,wave)
oldw<-readRDS(file.path(dd,"causal_participant_weights_internal_restricted.rds"))

# Over-identified CBPS is fit separately at each treatment decision. The
# denominator contains treatment history plus contemporaneous measured causes;
# the numerator contains treatment history only. WeightIt returns the ATE
# inverse generalized propensity weight, which is multiplied by the observed
# numerator probability to form a stabilized wave-specific weight.
denf<-fuel_state~previous_fuel+age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
denf0<-fuel_state~age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
obs_prob_num<-function(z){
 if(all(z$previous_fuel=="entry")){p<-prop.table(table(z$fuel_state));return(as.numeric(p[as.character(z$fuel_state)]))}
 fit<-multinom(fuel_state~previous_fuel,data=z,trace=FALSE,maxit=300)
 pr<-predict(fit,newdata=z,type="probs")
 if(is.null(dim(pr))){lev<-fit$lev;pr<-cbind(1-pr,pr);colnames(pr)<-lev}
 pr[cbind(seq_len(nrow(z)),match(as.character(z$fuel_state),colnames(pr)))]
}
res<-list(); audit<-list()
for(yr in sort(unique(long$wave))){
 cat("CBPS_START_WAVE=",yr,"\n",sep="");flush.console()
 z<-long%>%filter(wave==yr)%>%droplevels(); form<-if(all(z$previous_fuel=="entry"))denf0 else denf
 cb<-tryCatch(WeightIt::weightit(form,data=z,method="cbps",estimand="ATE",over=TRUE),error=function(e){cat("OVER_CBPS_FAILED_WAVE=",yr,"; ",conditionMessage(e),"\n",sep="");WeightIt::weightit(form,data=z,method="cbps",estimand="ATE",over=FALSE)})
 numerator<-obs_prob_num(z); z$cbps_wave_sw<-numerator*cb$weights
 z$cbps_wave_sw<-z$cbps_wave_sw/mean(z$cbps_wave_sw)
 if(any(!is.finite(z$cbps_wave_sw)))stop("Non-finite CBPS weight at ",yr)
 res[[as.character(yr)]]<-z
 audit[[as.character(yr)]]<-data.frame(wave=yr,n=nrow(z),states=n_distinct(z$fuel_state),minimum=min(z$cbps_wave_sw),P1=quantile(z$cbps_wave_sw,.01),P99=quantile(z$cbps_wave_sw,.99),maximum=max(z$cbps_wave_sw),ESS=sum(z$cbps_wave_sw)^2/sum(z$cbps_wave_sw^2),converged=TRUE)
 cat("CBPS_DONE_WAVE=",yr,"\n",sep="");flush.console()
}
cb_long<-bind_rows(res)%>%arrange(Idind,wave)%>%group_by(Idind)%>%mutate(cbps_treatment_weight_cum=cumprod(cbps_wave_sw))%>%ungroup()
cb_final<-cb_long%>%group_by(Idind)%>%summarise(cbps_treatment_weight=last(cbps_treatment_weight_cum),.groups="drop")%>%left_join(oldw%>%select(Idind,censoring_weight,strategy),by="Idind")%>%mutate(cbps_combined_weight=cbps_treatment_weight*censoring_weight)
saveRDS(cb_long,file.path(dd,"causal_CBPS_long_internal_restricted.rds"),compress="gzip")
saveRDS(cb_final,file.path(dd,"causal_CBPS_participant_weights_internal_restricted.rds"),compress="gzip")
write.csv(bind_rows(audit),file.path(od,"CBPS_wave_weight_audit.csv"),row.names=FALSE)

# Diagnostics use current three-state treatment at every wave. Pairwise maximum
# SMD is reported because a single binary SMD is invalid for multinomial fuel.
wmean<-function(x,w)sum(x*w,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)
wvar<-function(x,w){m<-wmean(x,w);sum(w*(x-m)^2,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)}
pair<-function(x,g,w,a,b){i<-g==a&is.finite(x);j<-g==b&is.finite(x);d<-sqrt((wvar(x[i],w[i])+wvar(x[j],w[j]))/2);if(!is.finite(d)||d==0)return(NA_real_);(wmean(x[j],w[j])-wmean(x[i],w[i]))/d}
maxsmd<-function(x,g,w){lev<-unique(as.character(g));ps<-combn(lev,2,simplify=FALSE);v<-vapply(ps,function(p)pair(x,as.character(g),w,p[1],p[2]),numeric(1));if(all(!is.finite(v)))NA_real_ else max(abs(v[is.finite(v)]))}
vars<-c("age","urban_i","education_i","income_i","assets_i","smoking_i","alcohol_i")
bal<-list()
for(yr in sort(unique(cb_long$wave))){z<-cb_long%>%filter(wave==yr);lo<-quantile(z$cbps_treatment_weight_cum,c(.01,.99));z$w99<-pmin(pmax(z$cbps_treatment_weight_cum,lo[1]),lo[2]);for(v in vars)bal[[length(bal)+1]]<-data.frame(wave=yr,variable=v,SMD_unweighted=maxsmd(z[[v]],z$fuel_state,rep(1,nrow(z))),SMD_CBPS=maxsmd(z[[v]],z$fuel_state,z$w99))}
balance<-bind_rows(bal)%>%mutate(abs_CBPS=abs(SMD_CBPS),meets_0_10=abs_CBPS<.10)
summary<-balance%>%group_by(wave)%>%summarise(max_abs_SMD=max(abs_CBPS,na.rm=TRUE),n_over_0_10=sum(abs_CBPS>=.10,na.rm=TRUE),.groups="drop")
write.csv(balance,file.path(td,"Table_causal_balance_CBPS.csv"),row.names=FALSE);write.csv(summary,file.path(od,"CBPS_balance_wave_summary.csv"),row.names=FALSE)

primary<-cb_final%>%filter(strategy%in%c("continued_any_solid","sustained_clean_only"));q<-quantile(primary$cbps_combined_weight,c(.01,.99));primary$cbps_combined_weight_99<-pmin(pmax(primary$cbps_combined_weight,q[1]),q[2])
ess<-primary%>%group_by(strategy)%>%summarise(N=n(),ESS=sum(cbps_combined_weight_99)^2/sum(cbps_combined_weight_99^2),.groups="drop")
write.csv(ess,file.path(od,"CBPS_primary_group_ESS.csv"),row.names=FALSE)
passed<-all(summary$max_abs_SMD<.10,na.rm=TRUE)&all(ess$ESS>=50)
 status<-data.frame(method="longitudinal CBPS-IPW (over-identified with exact-CBPS fallback for singular waves)",balance_pass=all(summary$max_abs_SMD<.10,na.rm=TRUE),ESS_pass=all(ess$ESS>=50),overall_pass=passed,decision=ifelse(passed,"Proceed to MSM with community-clustered inference","Do not fit causal outcome model; evaluate generalized overlap weighting with changed estimand"))
write.csv(status,file.path(od,"CBPS_final_status.csv"),row.names=FALSE)
print(summary);print(ess);print(status);cat("LONGITUDINAL_CBPS=",ifelse(passed,"PASS","FAIL"),"\n",sep="")
