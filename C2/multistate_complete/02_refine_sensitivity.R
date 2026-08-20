options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(splines);library(mgcv);library(sandwich)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);stopifnot(basename(project)=="longitudinal_upgrade")
d<-readRDS(file.path(project,"data","dynamic_analysis_internal_restricted.rds"));out<-file.path(project,"results")
rhs<-"ns(age,4)+sex+province_f+urban_f+education+income_asinh+income_missing+asset_count+asset_missing+smoking_f+alcohol_f"

# Correct GAM contrasts and pointwise SE using the lpmatrix covariance.
curves<-list();gd<-list()
for(y in c("HD_z","KDM_BAA")){
 f<-gam(as.formula(paste(y,"~s(solid_equiv_prop,k=4,bs='cr')+",rhs)),data=d,method="REML",na.action=na.exclude)
 nd<-d[rep(1,101),];nd$solid_equiv_prop<-seq(0,1,length.out=101)
 X<-predict(f,newdata=nd,type="lpmatrix");C<-X-matrix(X[1,],nrow(X),ncol(X),byrow=TRUE)
 est<-as.numeric(C%*%coef(f));se<-sqrt(pmax(rowSums((C%*%vcov(f))*C),0))
 curves[[y]]<-data.frame(outcome=y,solid_equiv_prop=nd$solid_equiv_prop,estimate=est,lower=est-1.96*se,upper=est+1.96*se)
 sm<-summary(f)$s.table;cv<-concurvity(f,full=TRUE);cv_est<-if(is.list(cv))cv$estimate else cv["estimate",]
 gd[[y]]<-data.frame(outcome=y,edf=sm[1,"edf"],F=sm[1,"F"],p_value=sm[1,"p-value"],deviance_explained=summary(f)$dev.expl,max_concurvity=max(cv_est[-1],na.rm=TRUE))
}
write.csv(bind_rows(curves),file.path(out,"gam_adjusted_curves.csv"),row.names=FALSE)
write.csv(bind_rows(gd),file.path(out,"gam_diagnostics.csv"),row.names=FALSE)

fit_one<-function(dat,outcome,expo,label){
 f<-lm(as.formula(paste(outcome,"~",expo,"+",rhs)),dat,na.action=na.exclude);used<-as.integer(rownames(model.frame(f)))
 V<-vcovCL(f,cluster=dat$commid[used],type="HC1",fix=TRUE);b<-coef(f)[expo];s<-sqrt(V[expo,expo])
 data.frame(analysis=label,outcome=outcome,term=expo,estimate=unname(b),lower=unname(b-1.96*s),upper=unname(b+1.96*s),p_value=unname(2*pnorm(-abs(b/s))),n=nobs(f),communities=n_distinct(dat$commid[used]))
}
sens<-bind_rows(
 fit_one(filter(d,n_intervals>=3),"HD_z","solid_prop10",">=3 observed intervals"),
 fit_one(filter(d,n_intervals>=3),"KDM_BAA","solid_prop10",">=3 observed intervals"),
 fit_one(filter(d,first_year<=1997),"HD_z","solid_prop10","history beginning by 1997"),
 fit_one(filter(d,first_year<=1997),"KDM_BAA","solid_prop10","history beginning by 1997"),
 fit_one(d,"HD_alt_z","solid_prop10","alternative HD construction"),
 fit_one(d,"KDM_BAA_direct_z","solid_prop10","direct KDM z-score"))
write.csv(sens,file.path(out,"dynamic_sensitivity_results.csv"),row.names=FALSE)
cat("REFINED_GAM_AND_SENSITIVITY=PASS\n")
