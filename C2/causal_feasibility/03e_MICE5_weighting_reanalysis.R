options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(mice);library(nnet);library(WeightIt);library(splines)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal");for(z in c(dd,od,td,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"03e_MICE5_weighting_reanalysis.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
base<-readRDS(file.path(dd,"causal_treatment_long_internal_restricted.rds"))%>%arrange(Idind,wave)
oldw<-readRDS(file.path(dd,"causal_participant_weights_internal_restricted.rds"))
waves<-sort(unique(base$wave));m<-5

# Wave-specific MICE avoids using future survey information to impute an earlier
# treatment-confounder history. Outcomes and biomarker variables are excluded.
completed<-replicate(m,base,simplify=FALSE);mi_audit<-list()
for(yr in waves){
 idx<-which(base$wave==yr);z<-base[idx,]
 impdat<-z%>%transmute(age,sex=factor(sex),province=factor(province),fuel_state=factor(fuel_state),previous_fuel=factor(previous_fuel),baseline_age_i,baseline_urban_i,baseline_education_i,baseline_income_i,baseline_assets_i,urban=factor(urban),education_years,income_asinh,asset_count,smoking=factor(smoking),alcohol=factor(alcohol))
 meth<-mice::make.method(impdat);meth[]<-"";targets<-c("urban","education_years","income_asinh","asset_count","smoking","alcohol")
 for(v in targets){if(anyNA(impdat[[v]])&&sum(!is.na(impdat[[v]]))>=50&&length(unique(impdat[[v]][!is.na(impdat[[v]])]))>1)meth[v]<-if(v%in%c("urban","smoking","alcohol"))"logreg" else "pmm"}
 pred<-mice::make.predictorMatrix(impdat);pred[,]<-0
 predictors<-setdiff(names(impdat),targets)
 for(v in targets)if(meth[v]!="")pred[v,setdiff(c(predictors,setdiff(targets,v)),v)]<-1
 diag(pred)<-0
 set.seed(20260815+yr);mids<-mice(impdat,m=m,maxit=5,method=meth,predictorMatrix=pred,printFlag=FALSE,seed=20260815+yr)
 for(j in seq_len(m)){
  cz<-complete(mids,j)
  # Entire-wave structural missingness remains represented by an indicator and
  # a neutral within-wave mode, exactly as in the complete-data pipeline.
  mode1<-function(x){u<-x[!is.na(x)];if(!length(u))return(0);as.numeric(names(sort(table(u),decreasing=TRUE)[1]))}
  to_num<-function(x)suppressWarnings(as.numeric(as.character(x)))
  completed[[j]]$urban_i[idx]<-ifelse(is.na(cz$urban),mode1(to_num(cz$urban)),to_num(cz$urban))
  completed[[j]]$education_i[idx]<-ifelse(is.na(cz$education_years),median(cz$education_years,na.rm=TRUE),cz$education_years)
  completed[[j]]$income_i[idx]<-ifelse(is.na(cz$income_asinh),median(cz$income_asinh,na.rm=TRUE),cz$income_asinh)
  completed[[j]]$assets_i[idx]<-ifelse(is.na(cz$asset_count),median(cz$asset_count,na.rm=TRUE),cz$asset_count)
  completed[[j]]$smoking_i[idx]<-ifelse(is.na(cz$smoking),mode1(to_num(cz$smoking)),to_num(cz$smoking))
  completed[[j]]$alcohol_i[idx]<-ifelse(is.na(cz$alcohol),mode1(to_num(cz$alcohol)),to_num(cz$alcohol))
 }
 mi_audit[[as.character(yr)]]<-data.frame(wave=yr,N=nrow(z),methods=paste(names(meth)[meth!=""],meth[meth!=""],sep="=",collapse=";"),structural_not_imputed=paste(names(meth)[names(meth)%in%targets&meth==""&vapply(impdat[names(meth)],anyNA,logical(1))],collapse=";"))
}
saveRDS(completed,file.path(dd,"causal_MICE5_long_internal_restricted.rds"),compress="gzip");write.csv(bind_rows(mi_audit),file.path(od,"MICE5_audit.csv"),row.names=FALSE)

num_prob<-function(z){if(all(z$previous_fuel=="entry")){p<-prop.table(table(z$fuel_state));return(as.numeric(p[as.character(z$fuel_state)]))};fit<-multinom(fuel_state~previous_fuel,data=z,trace=FALSE,maxit=300);pr<-predict(fit,newdata=z,type="probs");if(is.null(dim(pr))){pr<-cbind(1-pr,pr);colnames(pr)<-fit$lev};pr[cbind(seq_len(nrow(z)),match(as.character(z$fuel_state),colnames(pr)))]}
denf<-fuel_state~previous_fuel+ns(age,3)+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
denf0<-fuel_state~ns(age,3)+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
obsprob<-function(fit,z){pr<-predict(fit,newdata=z,type="probs");if(is.null(dim(pr))){pr<-cbind(1-pr,pr);colnames(pr)<-fit$lev};pr[cbind(seq_len(nrow(z)),match(as.character(z$fuel_state),colnames(pr))) ]}
wmean<-function(x,w)sum(x*w,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE);wvar<-function(x,w){q<-wmean(x,w);sum(w*(x-q)^2,na.rm=TRUE)/sum(w[is.finite(x)],na.rm=TRUE)}
pair<-function(x,g,w,a,b){i<-g==a&is.finite(x);j<-g==b&is.finite(x);d<-sqrt((wvar(x[i],w[i])+wvar(x[j],w[j]))/2);if(!is.finite(d)||d==0)return(NA_real_);(wmean(x[j],w[j])-wmean(x[i],w[i]))/d};maxsmd<-function(x,g,w){ps<-combn(unique(as.character(g)),2,simplify=FALSE);v<-vapply(ps,function(p)pair(x,as.character(g),w,p[1],p[2]),numeric(1));if(all(!is.finite(v)))NA_real_ else max(abs(v[is.finite(v)]))}
vars<-c("age","urban_i","education_i","income_i","assets_i","smoking_i","alcohol_i")
allbal<-list();alless<-list();weight_objects<-list()
for(j in seq_len(m)){
 d<-completed[[j]];parts<-list()
 for(yr in waves){z<-d%>%filter(wave==yr)%>%droplevels();f<-if(all(z$previous_fuel=="entry"))denf0 else denf
  fit<-multinom(f,data=z,trace=FALSE,maxit=500);z$mult_wave<-num_prob(z)/pmax(obsprob(fit,z),.001);z$mult_wave<-z$mult_wave/mean(z$mult_wave)
  ow<-WeightIt::weightit(f,data=z,method="glm",estimand="ATO");z$overlap_wave<-ow$weights/mean(ow$weights);parts[[as.character(yr)]]<-z}
 dl<-bind_rows(parts)%>%arrange(Idind,wave)%>%group_by(Idind)%>%mutate(mult_cum=cumprod(mult_wave),overlap_cum=cumprod(overlap_wave))%>%ungroup();weight_objects[[j]]<-dl
 for(yr in waves){z<-dl%>%filter(wave==yr);for(method in c("mult_cum","overlap_cum")){lo<-quantile(z[[method]],c(.01,.99));ww<-pmin(pmax(z[[method]],lo[1]),lo[2]);for(v in vars)allbal[[length(allbal)+1]]<-data.frame(imputation=j,wave=yr,method=method,variable=v,abs_SMD=maxsmd(z[[v]],z$fuel_state,ww))}}
 fin<-dl%>%group_by(Idind)%>%summarise(mult=last(mult_cum),overlap=last(overlap_cum),.groups="drop")%>%left_join(oldw%>%select(Idind,censoring_weight,strategy),by="Idind")%>%filter(strategy%in%c("continued_any_solid","sustained_clean_only"))
 for(method in c("mult","overlap")){x<-fin[[method]]*fin$censoring_weight;q<-quantile(x,c(.01,.99));x<-pmin(pmax(x,q[1]),q[2]);tmp<-fin%>%mutate(w=x)%>%group_by(strategy)%>%summarise(N=n(),ESS=sum(w)^2/sum(w^2),.groups="drop")%>%mutate(imputation=j,method=paste0(method,"_cum"));alless[[length(alless)+1]]<-tmp}
 cat("MICE_WEIGHTING_DONE=",j,"/5\n",sep="");flush.console()
}
saveRDS(weight_objects,file.path(dd,"causal_MICE5_weight_objects_internal_restricted.rds"),compress="gzip")
bal<-bind_rows(allbal);ess<-bind_rows(alless);summary<-bal%>%group_by(imputation,method,wave)%>%summarise(max_abs_SMD=max(abs_SMD,na.rm=TRUE),n_over_0_10=sum(abs_SMD>=.10,na.rm=TRUE),.groups="drop")
write.csv(bal,file.path(td,"Table_causal_balance_MICE5.csv"),row.names=FALSE);write.csv(summary,file.path(od,"MICE5_balance_wave_summary.csv"),row.names=FALSE);write.csv(ess,file.path(od,"MICE5_group_ESS.csv"),row.names=FALSE)
decision<-summary%>%group_by(method)%>%summarise(all_5_balance_pass=all(max_abs_SMD<.10),worst_max_abs_SMD=max(max_abs_SMD),.groups="drop")%>%left_join(ess%>%group_by(method)%>%summarise(all_5_ESS_pass=all(ESS>=50),minimum_ESS=min(ESS),.groups="drop"),by="method")%>%mutate(overall_pass=all_5_balance_pass&all_5_ESS_pass,MSM=ifelse(overall_pass,"permitted","withheld"))
write.csv(decision,file.path(od,"MICE5_weighting_final_status.csv"),row.names=FALSE);print(decision);cat("MICE5_REANALYSIS=COMPLETE\n")
