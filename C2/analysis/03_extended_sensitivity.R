options(stringsAsFactors=FALSE,cli.unicode=FALSE)
suppressPackageStartupMessages({library(dplyr);library(splines);library(sandwich)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="C2")stop("Run from C2")
out<-file.path(project,"results","formal")
d<-readRDS(file.path(project,"data","analysis_dataset_internal_restricted.rds"))
sel<-readRDS(file.path(project,"data","selection_cohort_internal_restricted.rds"))

vc<-function(fit,dd,cl="community")sandwich::vcovCL(fit,cluster=if(cl=="community")dd$commid else dd$hhid,type="HC1",fix=TRUE)
std<-function(outcome,rhs,dd,name){
 fit<-lm(as.formula(paste(outcome,"~",rhs)),data=dd);stopifnot(nobs(fit)==nrow(dd));V<-vc(fit,dd);n1<-dd;n0<-dd
 n1$trajectory<-factor("solid_to_clean_only",levels=levels(dd$trajectory));n0$trajectory<-factor("persistent_any_solid",levels=levels(dd$trajectory))
 tt<-delete.response(terms(fit));X1<-model.matrix(tt,n1,contrasts.arg=fit$contrasts,xlev=fit$xlevels);X0<-model.matrix(tt,n0,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
 cc<-colMeans(X1-X0);bb<-coef(fit);keep<-intersect(intersect(names(cc),names(bb)[!is.na(bb)]),colnames(V));e<-sum(cc[keep]*bb[keep]);s<-sqrt(as.numeric(t(cc[keep])%*%V[keep,keep,drop=FALSE]%*%cc[keep]))
 data.frame(outcome=outcome,analysis=name,estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_value=2*pnorm(-abs(e/s)),n=nrow(dd),clusters=n_distinct(dd$commid),cluster="community")
}
full_rhs<-"trajectory + ns(age,4) + sex + province_f + urban_f + education + income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
nosex_rhs<-"trajectory + ns(age,4) + province_f + urban_f + education + income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
p<-d%>%filter(trajectory%in%c("persistent_any_solid","solid_to_clean_only"))%>%droplevels()

# Frozen trajectory sensitivities.
p3<-p%>%filter(n_prior_effective>=3)%>%droplevels()
firstlast<-d%>%mutate(trajectory_fl=case_when(first_effective_state=="any_solid"&last_effective_state=="clean_only"~"solid_to_clean_only",
 first_effective_state=="any_solid"&last_effective_state=="any_solid"~"persistent_any_solid",TRUE~NA_character_))%>%filter(!is.na(trajectory_fl))
firstlast$trajectory<-factor(firstlast$trajectory_fl,levels=c("persistent_any_solid","solid_to_clean_only"))
sens<-bind_rows(std("HD_z",full_rhs,p3,">=3 prior effective fuel waves"),std("KDM_BAA",full_rhs,p3,">=3 prior effective fuel waves"),
 std("HD_z",full_rhs,firstlast,"First-last state definition"),std("KDM_BAA",full_rhs,firstlast,"First-last state definition"))

# Exploratory sex-stratified estimates and formal interaction tests.
for(sx in levels(p$sex)){
 z<-p%>%filter(sex==sx)%>%droplevels()
 sens<-bind_rows(sens,std("HD_z",nosex_rhs,z,paste0("Sex-stratified exploratory: ",sx)),std("KDM_BAA",nosex_rhs,z,paste0("Sex-stratified exploratory: ",sx)))
}
interaction_test<-function(outcome){
 rhs<-"trajectory*sex + ns(age,4) + province_f + urban_f + education + income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
 fit<-lm(as.formula(paste(outcome,"~",rhs)),data=p);V<-vc(fit,p);nm<-grep("trajectory.*:sex|sex.*:trajectory",names(coef(fit)),value=TRUE)
 stopifnot(length(nm)==1);e<-coef(fit)[nm];s<-sqrt(V[nm,nm]);data.frame(outcome=outcome,term=nm,interaction_estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_interaction=2*pnorm(-abs(e/s)),n=nrow(p),label="Exploratory sex interaction")
}
interactions<-bind_rows(interaction_test("HD_z"),interaction_test("KDM_BAA"));write.csv(interactions,file.path(out,"sex_interaction_tests.csv"),row.names=FALSE)

# Selection audit: formal inclusion among adults with all nine outcome inputs.
w_smd<-function(x,a,w=NULL){ok<-!is.na(x)&!is.na(a);x<-x[ok];a<-a[ok];if(is.null(w))w<-rep(1,length(x))else w<-w[ok];m1<-weighted.mean(x[a==1],w[a==1]);m0<-weighted.mean(x[a==0],w[a==0]);v1<-weighted.mean((x[a==1]-m1)^2,w[a==1]);v0<-weighted.mean((x[a==0]-m0)^2,w[a==0]);(m1-m0)/sqrt((v1+v0)/2)}
sel$selected_num<-as.integer(sel$selected_formal)
selbal<-bind_rows(lapply(c("age","urban","income_asinh","asset_count","smoking","alcohol","HD_z","KDM_BAA"),function(v)data.frame(variable=v,level="continuous/binary",smd_selected_vs_excluded=w_smd(sel[[v]],sel$selected_num))))
for(v in c("sex","education","province_f"))for(lv in levels(factor(sel[[v]])))selbal<-bind_rows(selbal,data.frame(variable=v,level=lv,smd_selected_vs_excluded=w_smd(as.numeric(sel[[v]]==lv),sel$selected_num)))
write.csv(selbal,file.path(out,"selection_balance.csv"),row.names=FALSE)
spfit<-glm(selected_num~ns(age,4)+sex+province_f+urban_f+education+income_asinh+income_missing+asset_count+asset_missing+smoking_f+alcohol_f,family=binomial(),data=sel)
stopifnot(nobs(spfit)==nrow(sel));sel$selection_ps<-pmin(pmax(predict(spfit,type="response"),.01),.99)
psmap<-sel%>%filter(selected_formal)%>%select(Idind,selection_ps);psens<-p%>%left_join(psmap,by="Idind")
psens$selection_weight<-1/psens$selection_ps;q<-quantile(psens$selection_weight,c(.01,.99));psens$selection_weight_trim<-pmin(pmax(psens$selection_weight,q[1]),q[2])
ess<-function(w)sum(w)^2/sum(w^2)
swd<-data.frame(n=nrow(psens),ps_min=min(psens$selection_ps),ps_p05=quantile(psens$selection_ps,.05),ps_median=median(psens$selection_ps),ps_p95=quantile(psens$selection_ps,.95),ps_max=max(psens$selection_ps),
 weight_min=min(psens$selection_weight),weight_sd=sd(psens$selection_weight),weight_p05=quantile(psens$selection_weight,.05),weight_p95=quantile(psens$selection_weight,.95),weight_p99=quantile(psens$selection_weight,.99),weight_max=max(psens$selection_weight),trim_lower=q[1],trim_upper=q[2],ESS_untrimmed=ess(psens$selection_weight),ESS_trimmed=ess(psens$selection_weight_trim),truncation="1st/99th percentile for exploratory selection-weight sensitivity")
write.csv(swd,file.path(out,"selection_weight_diagnostics.csv"),row.names=FALSE)
weighted_std<-function(outcome){fit<-lm(as.formula(paste(outcome,"~",full_rhs)),data=psens,weights=selection_weight_trim);V<-vc(fit,psens);n1<-psens;n0<-psens;n1$trajectory<-factor("solid_to_clean_only",levels=levels(psens$trajectory));n0$trajectory<-factor("persistent_any_solid",levels=levels(psens$trajectory));tt<-delete.response(terms(fit));X1<-model.matrix(tt,n1,contrasts.arg=fit$contrasts,xlev=fit$xlevels);X0<-model.matrix(tt,n0,contrasts.arg=fit$contrasts,xlev=fit$xlevels);cc<-colSums((X1-X0)*psens$selection_weight_trim)/sum(psens$selection_weight_trim);bb<-coef(fit);keep<-intersect(intersect(names(cc),names(bb)[!is.na(bb)]),colnames(V));e<-sum(cc[keep]*bb[keep]);s<-sqrt(as.numeric(t(cc[keep])%*%V[keep,keep,drop=FALSE]%*%cc[keep]));data.frame(outcome=outcome,analysis="Exploratory inverse selection weighted",estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_value=2*pnorm(-abs(e/s)),n=nrow(psens),clusters=n_distinct(psens$commid),cluster="community")}
sens<-bind_rows(sens,weighted_std("HD_z"),weighted_std("KDM_BAA"));write.csv(sens,file.path(out,"extended_sensitivity_results.csv"),row.names=FALSE)

# Individual deidentified reclassification output and aggregate matrix.
recl<-d%>%arrange(Idind)%>%transmute(record_id=sprintf("C2-%05d",row_number()),joint_primary_secondary=as.character(trajectory),only_primary=trajectory_primary_only,stacking_ever,reclassified=joint_primary_secondary!=only_primary)
write.csv(recl,file.path(out,"fuel_reclassification_individual_deidentified.csv"),row.names=FALSE)
recmat<-recl%>%count(only_primary,joint_primary_secondary,name="n")%>%group_by(only_primary)%>%mutate(row_percent=100*n/sum(n))%>%ungroup()
write.csv(recmat,file.path(out,"fuel_reclassification_matrix_counts.csv"),row.names=FALSE)

# A declared SAP deviation, not a silent omission.
dev<-data.frame(item="Multiple imputation of covariates",planned="Assess MI and compare with complete cases",implemented="Not run; primary retained missing categories/indicators and complete-case sensitivity was run",reason="The decision and its numerical diagnostics are recorded in the restricted internal analysis record.",impact="Recorded deviation; no claim that missing-category coding is equivalent to MI",status="Documented SAP deviation")
write.csv(dev,file.path(out,"sap_deviations.csv"),row.names=FALSE)
cat("PASS: extended trajectory, sex, selection and reclassification analyses written\n")
