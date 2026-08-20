options(stringsAsFactors=FALSE,warn=1,cli.unicode=FALSE)
suppressPackageStartupMessages({
  library(haven);library(dplyr);library(tidyr);library(splines);library(msm)
  library(sandwich);library(mgcv)
})
set.seed(20260813)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
if(basename(project)!="longitudinal_upgrade")stop("Run from C2/longitudinal_upgrade")
c2<-dirname(project); workspace<-dirname(c2)
out<-file.path(project,"results");logdir<-file.path(project,"logs")
dir.create(out,FALSE,TRUE);dir.create(logdir,FALSE,TRUE);dir.create(file.path(project,"data"),FALSE,TRUE)
logcon<-file(file.path(logdir,"01_multistate_dynamic_analysis.log"),"wt",encoding="UTF-8")
sink(logcon,split=TRUE);sink(logcon,type="message")
on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
cat("C2 longitudinal upgrade started",format(Sys.time()),"\n",R.version.string,"\n")

find_one<-function(pattern){
 z<-list.files(workspace,pattern=pattern,recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
 z<-z[grepl("CHNS-STATA",z,fixed=TRUE)]
 if(length(z)!=1)stop("Expected one file for ",pattern,"; found ",length(z))
 normalizePath(z,winslash="/",mustWork=TRUE)
}
asset_path<-find_one("^asset_12[.]dta$"); survey_path<-find_one("^surveys_pub_12[.]dta$")
asset<-read_dta(asset_path,col_select=c("hhid","wave","L8_1","L8_2"))%>%filter(wave<=2009)
survey<-read_dta(survey_path,col_select=c("Idind","hhid","wave"))%>%filter(wave<=2009)
ana<-readRDS(file.path(c2,"data","analysis_dataset_internal_restricted.rds"))

fuel_class<-function(x)case_when(x%in%c(1,6,7)~"solid",x%in%c(2,4,5)~"clean",TRUE~NA_character_)
panel<-survey%>%left_join(asset,by=c("hhid","wave"))%>%
 mutate(p=fuel_class(L8_1),s=fuel_class(L8_2),
  state_label=case_when(
   (p=="solid"|s=="solid") & !(p=="clean"|s=="clean")~"solid_only",
   (p=="solid"|s=="solid") & (p=="clean"|s=="clean")~"mixed",
   (p=="clean"|s=="clean") & !(p=="solid"|s=="solid")~"clean_only",
   TRUE~NA_character_),
  state=match(state_label,c("solid_only","mixed","clean_only")))%>%
 select(Idind,wave,state,state_label)%>%distinct()
stopifnot(!anyDuplicated(panel[c("Idind","wave")]))

# restrict multistate model to biomarker-linked adults, but use all their <=2009 observations
eligible_ids<-ana%>%pull(Idind)
mp<-panel%>%filter(Idind%in%eligible_ids,!is.na(state))%>%group_by(Idind)%>%
 filter(n()>=2,any(wave==2009))%>%arrange(wave,.by_group=TRUE)%>%ungroup()%>%
 left_join(ana%>%select(Idind,sex,urban_f)%>%distinct(),by="Idind")%>%
 mutate(time=wave-min(wave),urban_num=ifelse(as.character(urban_f)%in%c("Urban","1"),1,0),
        sex_num=ifelse(as.character(sex)%in%c("Male","1"),1,0))
if(n_distinct(mp$Idind)<500)stop("Too few multistate participants")
write.csv(mp%>%count(wave,state_label,name="n"),file.path(out,"state_counts_by_wave.csv"),row.names=FALSE)
flow<-data.frame(stage=c("2009 biomarker analysis cohort","At least two valid fuel states and valid 2009 state","Multistate observations"),
 n=c(length(unique(eligible_ids)),n_distinct(mp$Idind),nrow(mp)))
write.csv(flow,file.path(out,"sample_flow.csv"),row.names=FALSE)

# observed transitions and empirical rate initialization
pairs<-mp%>%group_by(Idind)%>%arrange(wave,.by_group=TRUE)%>%mutate(to=lead(state),to_label=lead(state_label),dt=lead(wave)-wave)%>%
 filter(!is.na(to),dt>0)%>%ungroup()
trans_counts<-pairs%>%count(state_label,to_label,name="n")
write.csv(trans_counts,file.path(out,"observed_transition_counts.csv"),row.names=FALSE)
ptime<-pairs%>%group_by(state)%>%summarise(py=sum(dt),.groups="drop")
Q<-matrix(.01,3,3);diag(Q)<-0
for(i in 1:3)for(j in 1:3)if(i!=j){nn<-sum(pairs$state==i & pairs$to==j);py<-ptime$py[ptime$state==i];Q[i,j]<-max(nn/py,.002)}
diag(Q)<--rowSums(Q)

fit0<-msm(state~wave,subject=Idind,data=mp,qmatrix=Q,center=FALSE,control=list(fnscale=4000,maxit=20000))
if(!isTRUE(fit0$opt$convergence==0))stop("Base msm did not converge")
qe<-qmatrix.msm(fit0,ci="normal")
qtab<-data.frame(from=rep(c("solid_only","mixed","clean_only"),each=3),to=rep(c("solid_only","mixed","clean_only"),3),
 estimate=as.vector(t(qe$estimates)),lower=as.vector(t(qe$L)),upper=as.vector(t(qe$U)))%>%filter(from!=to)
write.csv(qtab,file.path(out,"transition_intensities.csv"),row.names=FALSE)

p5<-pmatrix.msm(fit0,t=5);p10<-pmatrix.msm(fit0,t=10)
pmat_df<-function(x,t){z<-as.data.frame(unclass(x));names(z)<-c("solid_only","mixed","clean_only");z%>%mutate(from=c("solid_only","mixed","clean_only"))%>%
 pivot_longer(-from,names_to="to",values_to="estimate")%>%mutate(years=t)%>%select(years,from,to,estimate)}
pmats<-bind_rows(pmat_df(p5,5),pmat_df(p10,10));write.csv(pmats,file.path(out,"transition_probabilities.csv"),row.names=FALSE)
soj<-sojourn.msm(fit0);sojdf<-data.frame(state=c("solid_only","mixed","clean_only"),estimate=soj$estimates,lower=soj$L,upper=soj$U)
write.csv(sojdf,file.path(out,"mean_sojourn_years.csv"),row.names=FALSE)

occ<-bind_rows(lapply(seq(0,20,.25),function(tt){P<-pmatrix.msm(fit0,t=tt);data.frame(year=tt,state=c("solid_only","mixed","clean_only"),probability=P[1,])}))
write.csv(occ,file.path(out,"occupancy_from_solid.csv"),row.names=FALSE)

# exploratory common covariate effects across transitions
fitcov<-suppressWarnings(try(msm(state~wave,subject=Idind,data=mp,qmatrix=Q,covariates=~sex_num+urban_num,center=FALSE,control=list(fnscale=4000,maxit=20000)),silent=TRUE))
cov_ok<-!inherits(fitcov,"try-error")&&fitcov$opt$convergence==0&&all(eigen(fitcov$opt$hessian,symmetric=TRUE,only.values=TRUE)$values>0)
if(cov_ok){hr<-hazard.msm(fitcov);capture.output(hr,file=file.path(out,"transition_covariate_HR.txt"))
}else writeLines("Exploratory covariate msm had an unstable/non-positive-definite Hessian; estimates were not reported and the base Q model was retained.",file.path(out,"transition_covariate_HR.txt"))

# individual interval histories, left-state allocation
histories<-panel%>%filter(Idind%in%eligible_ids,!is.na(state))%>%group_by(Idind)%>%arrange(wave,.by_group=TRUE)%>%
 filter(n()>=2,any(wave==2009))%>%mutate(next_wave=lead(wave),next_state=lead(state),dt=next_wave-wave)%>%
 filter(!is.na(dt),dt>0,next_wave<=2009)%>%ungroup()
summ<-histories%>%group_by(Idind)%>%summarise(
 observed_years=sum(dt),solid_years=sum(dt*(state==1)),mixed_years=sum(dt*(state==2)),clean_years=sum(dt*(state==3)),
 solid_equiv_years=solid_years+.5*mixed_years,solid_equiv_prop=solid_equiv_years/observed_years,
 n_intervals=n(),n_changes=sum(state!=next_state),n_reversals=sum(state==3 & next_state%in%c(1,2)),
 any_reversal=as.integer(n_reversals>0),first_year=min(wave),.groups="drop")
# sustained clean duration: last continuous run ending in clean at 2009, based on observed intervals
cleanrun<-panel%>%filter(Idind%in%eligible_ids,!is.na(state))%>%group_by(Idind)%>%arrange(wave,.by_group=TRUE)%>%
 filter(n()>=2,any(wave==2009))%>%summarise(
 sustained_clean_years={w<-wave;s<-state;if(tail(w,1)!=2009||tail(s,1)!=3)0 else {k<-length(s);while(k>1&&s[k-1]==3)k<-k-1;2009-w[k]}},.groups="drop")
summ<-left_join(summ,cleanrun,by="Idind")
stopifnot(all(abs(summ$solid_years+summ$mixed_years+summ$clean_years-summ$observed_years)<1e-8),all(summ$solid_equiv_prop>=0&summ$solid_equiv_prop<=1))

d<-ana%>%left_join(summ,by="Idind")%>%filter(!is.na(solid_equiv_prop))
d<-d%>%mutate(solid_prop10=solid_equiv_prop/.10,clean5=sustained_clean_years/5,rev_any=any_reversal)
saveRDS(d,file.path(project,"data","dynamic_analysis_internal_restricted.rds"),compress="gzip")
deid<-d%>%arrange(Idind)%>%mutate(record_id=sprintf("C2L-%05d",row_number()))%>%select(record_id,everything(),-any_of(c("Idind","hhid","commid","stratum","Idind_key")))
saveRDS(deid,file.path(project,"data","dynamic_analysis_deidentified.rds"),compress="gzip")
write.csv(summ%>%summarise(n=n(),median_observed_years=median(observed_years),p25=q25<-quantile(observed_years,.25),p75=quantile(observed_years,.75),
 mean_solid_equiv_prop=mean(solid_equiv_prop),reversal_n=sum(any_reversal),sustained_clean_n=sum(sustained_clean_years>0)),file.path(out,"dynamic_history_summary.csv"),row.names=FALSE)
write.csv(d%>%count(cut(solid_equiv_prop,breaks=seq(0,1,.1),include.lowest=TRUE),name="n"),file.path(out,"exposure_support_deciles.csv"),row.names=FALSE)

rhs<-"ns(age,4)+sex+province_f+urban_f+education+income_asinh+income_missing+asset_count+asset_missing+smoking_f+alcohol_f"
fit_extract<-function(outcome,expo,label,unit,dat=d){
 f<-lm(as.formula(paste(outcome,"~",expo,"+",rhs)),data=dat,na.action=na.exclude)
 used<-as.integer(rownames(model.frame(f)));V<-sandwich::vcovCL(f,cluster=dat$commid[used],type="HC1",fix=TRUE)
 b<-coef(f)[expo];se<-sqrt(V[expo,expo]);data.frame(outcome=outcome,exposure=label,unit=unit,estimate=unname(b),lower=unname(b-1.96*se),upper=unname(b+1.96*se),p_value=unname(2*pnorm(-abs(b/se))),n=nobs(f),communities=n_distinct(dat$commid[used]))
}
res<-bind_rows(
 fit_extract("HD_z","solid_prop10","Solid-equivalent history","per 10 percentage points"),
 fit_extract("KDM_BAA","solid_prop10","Solid-equivalent history","per 10 percentage points"),
 fit_extract("HD_z","clean5","Sustained clean-only duration","per 5 years"),
 fit_extract("KDM_BAA","clean5","Sustained clean-only duration","per 5 years"),
 fit_extract("HD_z","rev_any","Clean-to-solid/mixed reversal","any versus none"),
 fit_extract("KDM_BAA","rev_any","Clean-to-solid/mixed reversal","any versus none"))
write.csv(res,file.path(out,"dynamic_association_results.csv"),row.names=FALSE)

# GAM nonlinear sensitivity, using ML for comparable EDF/P values
gamfits<-list();gamdiag<-list();curves<-list()
for(y in c("HD_z","KDM_BAA")){
 f<-gam(as.formula(paste(y,"~s(solid_equiv_prop,k=4,bs='cr')+",rhs)),data=d,method="REML",na.action=na.exclude)
 gamfits[[y]]<-f;sm<-summary(f)$s.table
 cv<-concurvity(f,full=TRUE);cv_est<-if(is.list(cv))cv$estimate else cv["estimate",]
 gamdiag[[y]]<-data.frame(outcome=y,edf=sm[1,"edf"],F=sm[1,"F"],p_value=sm[1,"p-value"],deviance_explained=summary(f)$dev.expl,max_concurvity=max(cv_est[-1],na.rm=TRUE))
 nd<-d[rep(1,101),];nd$solid_equiv_prop<-seq(0,1,length.out=101)
 pr<-predict(f,newdata=nd,se.fit=TRUE,type="response",exclude=NULL)
 # center curve at zero exposure for adjusted within-profile contrast
 base<-pr$fit[1]; curves[[y]]<-data.frame(outcome=y,solid_equiv_prop=nd$solid_equiv_prop,estimate=pr$fit-base,lower=pr$fit-base-1.96*pr$se.fit,upper=pr$fit-base+1.96*pr$se.fit)
}
write.csv(bind_rows(gamdiag),file.path(out,"gam_diagnostics.csv"),row.names=FALSE)
write.csv(bind_rows(curves),file.path(out,"gam_adjusted_curves.csv"),row.names=FALSE)

# community bootstrap of frozen-summary regression contrasts
communities<-unique(d$commid);B<-1000L;boot<-vector("list",B)
for(bi in seq_len(B)){
 sam<-sample(communities,length(communities),replace=TRUE)
 db<-bind_rows(lapply(seq_along(sam),function(k)d[d$commid==sam[k],]%>%mutate(comm_boot=k)))
 one<-lapply(list(c("HD_z","solid_prop10"),c("KDM_BAA","solid_prop10"),c("HD_z","clean5"),c("KDM_BAA","clean5"),c("HD_z","rev_any"),c("KDM_BAA","rev_any")),function(z){
   f<-lm(as.formula(paste(z[1],"~",z[2],"+",rhs)),data=db);data.frame(outcome=z[1],term=z[2],estimate=unname(coef(f)[z[2]]))})
 boot[[bi]]<-bind_rows(one)%>%mutate(replicate=bi)
 if(bi%%100==0)cat("bootstrap",bi,"/",B,"\n")
}
boot<-bind_rows(boot);write.csv(boot,file.path(out,"community_bootstrap_replicates.csv"),row.names=FALSE)
bs<-boot%>%group_by(outcome,term)%>%summarise(lower=quantile(estimate,.025,na.rm=TRUE,type=8),upper=quantile(estimate,.975,na.rm=TRUE,type=8),successful=sum(is.finite(estimate)),.groups="drop")
write.csv(bs,file.path(out,"community_bootstrap_summary.csv"),row.names=FALSE)

# simple diagnostics for primary models
diagout<-bind_rows(lapply(c("HD_z","KDM_BAA"),function(y){f<-lm(as.formula(paste(y,"~solid_prop10+",rhs)),d);data.frame(outcome=y,n=nobs(f),adj_r2=summary(f)$adj.r.squared,max_cooks=max(cooks.distance(f)),resid_skew=mean((residuals(f)-mean(residuals(f)))^3)/sd(residuals(f))^3)}))
write.csv(diagout,file.path(out,"regression_diagnostics.csv"),row.names=FALSE)
cat("MULTISTATE_DYNAMIC_ANALYSIS=PASS\n")
cat("MSM_N=",n_distinct(mp$Idind)," OUTCOME_N=",nrow(d)," BOOTSTRAP=",B,"\n",sep="")
