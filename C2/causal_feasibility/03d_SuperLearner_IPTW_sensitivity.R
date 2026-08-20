options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(SuperLearner);library(nnet)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal");for(z in c(dd,od,td,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"03d_SuperLearner_IPTW_sensitivity.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
long<-readRDS(file.path(dd,"causal_treatment_long_internal_restricted.rds"))%>%arrange(Idind,wave)
oldw<-readRDS(file.path(dd,"causal_participant_weights_internal_restricted.rds"))
learners<-c("SL.mean","SL.glm","SL.gam");if(requireNamespace("glmnet",quietly=TRUE))learners<-c(learners,"SL.glmnet");if(requireNamespace("ranger",quietly=TRUE))learners<-c(learners,"SL.ranger")
cat("SL_LIBRARY=",paste(learners,collapse=","),"\n",sep="")
numprob<-function(z){if(all(z$previous_fuel=="entry")){p<-prop.table(table(z$fuel_state));return(as.numeric(p[as.character(z$fuel_state)]))};fit<-multinom(fuel_state~previous_fuel,data=z,trace=FALSE,maxit=300);pr<-predict(fit,newdata=z,type="probs");if(is.null(dim(pr))){pr<-cbind(1-pr,pr);colnames(pr)<-fit$lev};pr[cbind(seq_len(nrow(z)),match(as.character(z$fuel_state),colnames(pr)))]}
out<-list();aud<-list()
for(yr in sort(unique(long$wave))){
 z<-long%>%filter(wave==yr)%>%droplevels();lev<-levels(z$fuel_state)
 xf<-if(all(z$previous_fuel=="entry"))~age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m else ~previous_fuel+age+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
 X<-model.matrix(xf,data=z)[,-1,drop=FALSE]
 keep<-apply(X,2,function(x)sd(x)>0);X<-as.data.frame(X[,keep,drop=FALSE]);pred<-matrix(NA_real_,nrow(z),length(lev),dimnames=list(NULL,lev));risks<-list()
 for(k in seq_along(lev)){set.seed(20260815+yr+k);sl<-SuperLearner(Y=as.integer(z$fuel_state==lev[k]),X=X,family=binomial(),SL.library=learners,cvControl=list(V=5),verbose=FALSE);pred[,k]<-sl$SL.predict;risks[[k]]<-data.frame(wave=yr,state=lev[k],learner=names(sl$coef),coefficient=as.numeric(sl$coef),cv_risk=as.numeric(sl$cvRisk))}
 pred<-pmax(pmin(pred,.999),.001);pred<-pred/rowSums(pred);den<-pred[cbind(seq_len(nrow(z)),match(as.character(z$fuel_state),lev))]
 z$sl_wave_sw<-numprob(z)/den;z$sl_wave_sw<-z$sl_wave_sw/mean(z$sl_wave_sw);out[[as.character(yr)]]<-z;aud[[as.character(yr)]]<-bind_rows(risks)
 cat("SL_DONE_WAVE=",yr,"\n",sep="");flush.console()
}
sllong<-bind_rows(out)%>%arrange(Idind,wave)%>%group_by(Idind)%>%mutate(sl_treatment_weight_cum=cumprod(sl_wave_sw))%>%ungroup()
slfinal<-sllong%>%group_by(Idind)%>%summarise(sl_treatment_weight=last(sl_treatment_weight_cum),.groups="drop")%>%left_join(oldw%>%select(Idind,censoring_weight,strategy),by="Idind")%>%mutate(sl_combined_weight=sl_treatment_weight*censoring_weight)
saveRDS(sllong,file.path(dd,"causal_SL_IPTW_long_internal_restricted.rds"),compress="gzip");saveRDS(slfinal,file.path(dd,"causal_SL_IPTW_participant_weights_internal_restricted.rds"),compress="gzip");write.csv(bind_rows(aud),file.path(od,"SL_learner_audit.csv"),row.names=FALSE)

wmean<-function(x,w)sum(x*w,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE);wvar<-function(x,w){m<-wmean(x,w);sum(w*(x-m)^2,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)}
pair<-function(x,g,w,a,b){i<-g==a&is.finite(x);j<-g==b&is.finite(x);d<-sqrt((wvar(x[i],w[i])+wvar(x[j],w[j]))/2);if(!is.finite(d)||d==0)return(NA_real_);(wmean(x[j],w[j])-wmean(x[i],w[i]))/d};maxsmd<-function(x,g,w){ps<-combn(unique(as.character(g)),2,simplify=FALSE);v<-vapply(ps,function(p)pair(x,as.character(g),w,p[1],p[2]),numeric(1));if(all(!is.finite(v)))NA_real_ else max(abs(v[is.finite(v)]))}
vars<-c("age","urban_i","education_i","income_i","assets_i","smoking_i","alcohol_i");bal<-list()
for(yr in sort(unique(sllong$wave))){z<-sllong%>%filter(wave==yr);lo<-quantile(z$sl_treatment_weight_cum,c(.01,.99));z$w99<-pmin(pmax(z$sl_treatment_weight_cum,lo[1]),lo[2]);for(v in vars)bal[[length(bal)+1]]<-data.frame(wave=yr,variable=v,SMD_SL=maxsmd(z[[v]],z$fuel_state,z$w99))}
balance<-bind_rows(bal)%>%mutate(abs_SL=abs(SMD_SL),meets_0_10=abs_SL<.10);bsum<-balance%>%group_by(wave)%>%summarise(max_abs_SMD=max(abs_SL,na.rm=TRUE),n_over_0_10=sum(abs_SL>=.10,na.rm=TRUE),.groups="drop")
write.csv(balance,file.path(td,"Table_causal_balance_SL_IPTW.csv"),row.names=FALSE);write.csv(bsum,file.path(od,"SL_balance_wave_summary.csv"),row.names=FALSE)
primary<-slfinal%>%filter(strategy%in%c("continued_any_solid","sustained_clean_only"));q<-quantile(primary$sl_combined_weight,c(.01,.99));primary$w99<-pmin(pmax(primary$sl_combined_weight,q[1]),q[2]);ess<-primary%>%group_by(strategy)%>%summarise(N=n(),ESS=sum(w99)^2/sum(w99^2),.groups="drop");write.csv(ess,file.path(od,"SL_primary_group_ESS.csv"),row.names=FALSE)
pass<-all(bsum$max_abs_SMD<.10,na.rm=TRUE)&&all(ess$ESS>=50);status<-data.frame(method="one-vs-rest Super Learner stabilized IPTW sensitivity",role="sensitivity only",balance_pass=all(bsum$max_abs_SMD<.10,na.rm=TRUE),ESS_pass=all(ess$ESS>=50),overall_pass=pass,decision=ifelse(pass,"May report as model-dependence sensitivity; not promoted automatically to primary","Does not rescue causal identification"));write.csv(status,file.path(od,"SL_IPTW_final_status.csv"),row.names=FALSE)
print(bsum);print(ess);print(status);cat("SUPER_LEARNER_IPTW_SENSITIVITY=",ifelse(pass,"PASS","FAIL"),"\n",sep="")
