options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(mice)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);dd<-file.path(project,"data","causal_landmark");od<-file.path(project,"results","causal_landmark");ld<-file.path(project,"logs","causal_landmark")
d<-readRDS(file.path(dd,"landmark_analysis_no_outcomes_internal_restricted.rds"))
impvars<-c("Idind","A","age2004","sex","province","fuel2004","previous_fuel_state","number_prior_fuel_waves","proportion_prior_solid_only","proportion_prior_mixed","ever_clean_before_2006","ever_mixed_before_2006","number_fuel_transitions","number_clean_to_mixed_reversals","number_mixed_to_solid_reversals","years_since_first_clean_observation","duration_weighted_solid_equivalent_history_to_2004","last_two_wave_pattern","urban_recent","education_recent","income_recent","income_mean","income_change","assets_recent","assets_mean","assets_change","smoking_recent","alcohol_recent","commid2004","biomarker2009","fuel2009")
x<-d[,impvars];x$A<-factor(x$A);x$urban_recent<-factor(x$urban_recent);x$smoking_recent<-factor(x$smoking_recent);x$alcohol_recent<-factor(x$alcohol_recent)
meth<-make.method(x);meth[]<-"";cont<-c("age2004","education_recent","income_recent","income_mean","income_change","assets_recent","assets_mean","assets_change","years_since_first_clean_observation");bin<-c("urban_recent","smoking_recent","alcohol_recent")
for(v in cont)if(anyNA(x[[v]]))meth[v]<-"pmm";for(v in bin)if(anyNA(x[[v]])&&length(unique(x[[v]][!is.na(x[[v]])]))>1)meth[v]<-"logreg"
pred<-make.predictorMatrix(x);pred[,]<-0;allowed<-setdiff(names(x),c("Idind","biomarker2009","fuel2009"));for(v in names(meth)[meth!=""])pred[v,setdiff(allowed,v)]<-1;diag(pred)<-0
set.seed(20260815);mi<-mice(x,m=5,maxit=5,method=meth,predictorMatrix=pred,printFlag=FALSE,seed=20260815)
completed<-lapply(1:5,function(i){z<-complete(mi,i);z$A<-as.integer(as.character(z$A));z})
saveRDS(completed,file.path(dd,"landmark_MICE5_no_outcomes_internal_restricted.rds"),compress="gzip")
write.csv(data.frame(variable=names(meth),method=meth,missingness=vapply(x,function(z)mean(is.na(z)),numeric(1))),file.path(od,"landmark_MICE5_audit.csv"),row.names=FALSE);writeLines("LANDMARK_MICE5=PASS_NO_OUTCOMES",file.path(ld,"03_landmark_MICE.log"),useBytes=TRUE)
