options(stringsAsFactors=FALSE, cli.unicode=FALSE)
suppressPackageStartupMessages({library(dplyr);library(splines);library(sandwich)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="C2")stop("Run from C2")
out<-file.path(project,"results","formal");d<-readRDS(file.path(project,"data","analysis_dataset.rds"))
p<-d%>%filter(trajectory%in%c("persistent_any_solid","solid_to_clean_only"))%>%droplevels();p$treated<-as.integer(p$trajectory=="solid_to_clean_only")
full_rhs<-"trajectory + ns(age,4) + sex + province_f + urban_f + education + income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
min_rhs<-"trajectory + ns(age,4) + sex + province_f + urban_f"
vc<-function(fit,dd,cl="community")sandwich::vcovCL(fit,cluster=if(cl=="community")dd$commid else dd$hhid,type="HC1",fix=TRUE)
std<-function(rhs,name,cl="community"){
 fit<-lm(as.formula(paste("KDM_BAA~",rhs)),data=p);stopifnot(nobs(fit)==nrow(p));V<-vc(fit,p,cl);n1<-p;n0<-p;n1$trajectory<-factor("solid_to_clean_only",levels=levels(p$trajectory));n0$trajectory<-factor("persistent_any_solid",levels=levels(p$trajectory));tt<-delete.response(terms(fit));X1<-model.matrix(tt,n1,contrasts.arg=fit$contrasts,xlev=fit$xlevels);X0<-model.matrix(tt,n0,contrasts.arg=fit$contrasts,xlev=fit$xlevels);cc<-colMeans(X1-X0);bb<-coef(fit);keep<-intersect(intersect(names(cc),names(bb)[!is.na(bb)]),colnames(V));e<-sum(cc[keep]*bb[keep]);s<-sqrt(as.numeric(t(cc[keep])%*%V[keep,keep,drop=FALSE]%*%cc[keep]));data.frame(outcome="KDM_BAA",analysis=name,contrast="solid_to_clean_only vs persistent_any_solid",estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_value=2*pnorm(-abs(e/s)),n=nrow(p),clusters=dplyr::n_distinct(if(cl=="community")p$commid else p$hhid),cluster=cl)
}
psfit<-glm(treated~ns(age,4)+sex+province_f+urban_f+education+income_asinh+income_missing+asset_count+asset_missing+smoking_f+alcohol_f,family=binomial(),data=p);stopifnot(nobs(psfit)==nrow(p));p$ps<-pmin(pmax(predict(psfit,type="response"),1e-4),1-1e-4);p$ow<-ifelse(p$treated==1,1-p$ps,p$ps)
wfit<-lm(KDM_BAA~treated,data=p,weights=ow);V<-vc(wfit,p);e<-coef(wfit)["treated"];s<-sqrt(V["treated","treated"])
wr<-data.frame(outcome="KDM_BAA",analysis="Overlap weighted",contrast="solid_to_clean_only vs persistent_any_solid",estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_value=2*pnorm(-abs(e/s)),n=nrow(p),clusters=n_distinct(p$commid),cluster="community")
yr<-bind_rows(std("trajectory","Unadjusted"),std(min_rhs,"Minimally adjusted"),std(full_rhs,"Fully adjusted"),std(full_rhs,"Fully adjusted; household clustered","household"),wr)
main<-read.csv(file.path(out,"main_secondary_results.csv"));main<-main[main$outcome!="KDM_BAA",];write.csv(bind_rows(main,yr),file.path(out,"main_secondary_results.csv"),row.names=FALSE)
raw_missing<-data.frame(variable=c("education_level_code","hhincpc_cpi","asset_count_raw","smoking","alcohol"),n_missing=c(sum(is.na(p$education_level_code)),sum(is.na(p$hhincpc_cpi)),sum(p$asset_missing==1),sum(is.na(p$smoking)),sum(is.na(p$alcohol))))
raw_missing$pct_missing=100*raw_missing$n_missing/nrow(p);write.csv(raw_missing,file.path(out,"raw_covariate_missingness_primary_contrast.csv"),row.names=FALSE)
integrity<-data.frame(check=c("nominal primary n","full outcome model n","propensity model n","KDM nonmissing"),value=c(nrow(p),nobs(lm(as.formula(paste("KDM_BAA~",full_rhs)),data=p)),nobs(psfit),sum(!is.na(p$KDM_BAA))),pass=c(TRUE,nobs(lm(as.formula(paste("KDM_BAA~",full_rhs)),data=p))==nrow(p),nobs(psfit)==nrow(p),sum(!is.na(p$KDM_BAA))==nrow(p)))
write.csv(integrity,file.path(out,"model_n_integrity.csv"),row.names=FALSE);cat("PASS: raw-year KDM and model-n integrity outputs written; n=",nrow(p),"\n",sep="")
