options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(nnet);library(splines)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");ld<-file.path(project,"logs","causal")
for(z in c(dd,od,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"03_IPTW_IPCW.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
d0<-readRDS(file.path(dd,"causal_long_full_internal_restricted.rds"))
support<-readRDS(file.path(dd,"target_strategy_support_internal_restricted.rds"))
waves<-c(1989,1991,1993,1997,2000,2004,2006,2009)

# Imputation is wave-specific median/mode plus explicit missing indicators; no future values are used.
mode1<-function(x){z<-x[!is.na(x)];if(!length(z))return(0);as.numeric(names(sort(table(z),decreasing=TRUE)[1]))}
prep<-function(d){d%>%group_by(wave)%>%mutate(
 urban_m=is.na(urban),education_m=is.na(education_years),income_m=is.na(income_asinh),assets_m=is.na(asset_count),smoking_m=is.na(smoking),alcohol_m=is.na(alcohol),
 urban_i=ifelse(is.na(urban),mode1(urban),urban),education_i=ifelse(is.na(education_years),median(education_years,na.rm=TRUE),education_years),
 income_i=ifelse(is.na(income_asinh),median(income_asinh,na.rm=TRUE),income_asinh),assets_i=ifelse(is.na(asset_count),median(asset_count,na.rm=TRUE),asset_count),
 smoking_i=ifelse(is.na(smoking),mode1(smoking),smoking),alcohol_i=ifelse(is.na(alcohol),mode1(alcohol),alcohol))%>%ungroup()%>%
 mutate(sex=factor(sex),province_f=factor(province),fuel_state=factor(fuel_state,levels=c("solid_only","mixed","clean_only")),
  previous_fuel=ifelse(is.na(previous_fuel),"entry",as.character(previous_fuel)),previous_fuel=factor(previous_fuel,levels=c("entry","solid_only","mixed","clean_only")))}
d<-prep(d0)
complete_ids<-support%>%filter(complete_sequence)%>%pull(Idind)

# Baseline missingness handled once, using cohort medians/modes.
for(v in c("baseline_age","baseline_urban","baseline_education","baseline_income","baseline_assets")){
 d[[paste0(v,"_m")]]<-is.na(d[[v]]); med<-median(d[[v]],na.rm=TRUE);if(!is.finite(med))med<-0;d[[paste0(v,"_i")]]<-ifelse(is.na(d[[v]]),med,d[[v]])
}

fit_state_model<-function(form,data){z<-droplevels(data);lev<-levels(z$fuel_state);if(length(lev)<2)stop("Only one treatment state in a wave");
 if(length(lev)==2)list(fit=glm(form,data=z,family=binomial()),lev=lev,type="binary") else list(fit=multinom(form,data=z,trace=FALSE,maxit=500),lev=lev,type="multi")}
get_obs_prob<-function(obj,newd,state){if(obj$type=="binary"){p<-predict(obj$fit,newdata=newd,type="response");pr<-cbind(1-p,p);colnames(pr)<-obj$lev}else{pr<-predict(obj$fit,newdata=newd,type="probs");if(is.null(dim(pr)))pr<-matrix(pr,ncol=length(obj$lev));colnames(pr)<-obj$lev};idx<-match(as.character(state),colnames(pr));out<-pr[cbind(seq_len(nrow(newd)),idx)];pmax(pmin(out,.999),.001)}
num_formula<-fuel_state~previous_fuel
# Prespecified flexible refinement: nonlinear socioeconomic terms and interactions
# with prior treatment capture the empirically strong path dependence without
# introducing future information or outcome-driven variable selection.
den_formula<-fuel_state~previous_fuel*(ns(age,3)+urban_i+ns(education_i,3)+ns(income_i,3)+ns(assets_i,3)+smoking_i+alcohol_i)+sex+province_f+province_f:urban_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
den_entry_formula<-fuel_state~ns(age,3)+urban_i+ns(education_i,3)+ns(income_i,3)+ns(assets_i,3)+smoking_i+alcohol_i+sex+province_f+province_f:urban_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
d$tp_num<-d$tp_den<-NA_real_
model_audit<-list()
for(w in waves){
 # The treatment model targets the linked complete-history trial cohort defined in
 # the frozen target-trial specification. Attrition/biomarker selection are
 # handled separately by IPCW below; mixing non-trial entrants into this model
 # changes the propensity-score population and degraded covariate balance.
 fitd<-d%>%filter(Idind%in%complete_ids,wave==w,!is.na(fuel_state));if(n_distinct(fitd$fuel_state)<2)stop("Only one treatment state at wave ",w)
 nf<-if(all(fitd$previous_fuel=="entry"))update(num_formula,.~.-previous_fuel)else num_formula
 df<-if(all(fitd$previous_fuel=="entry"))den_entry_formula else den_formula
 fn<-fit_state_model(nf,fitd);fd<-fit_state_model(df,fitd)
 ii<-which(d$Idind%in%complete_ids&d$wave==w&!is.na(d$fuel_state));d$tp_num[ii]<-get_obs_prob(fn,d[ii,],d$fuel_state[ii]);d$tp_den[ii]<-get_obs_prob(fd,d[ii,],d$fuel_state[ii])
 model_audit[[as.character(w)]]<-data.frame(wave=w,n=nrow(fitd),states=paste(names(table(fitd$fuel_state)),table(fitd$fuel_state),collapse=";"),min_den=min(d$tp_den[ii]),max_den=max(d$tp_den[ii]))
}
write.csv(bind_rows(model_audit),file.path(od,"treatment_model_audit.csv"),row.names=FALSE)

# Complete observed treatment histories reaching 2009 define the primary per-protocol cohort.
outcome_rows<-d%>%filter(Idind%in%complete_ids,biomarker2009==1,!is.na(fuel_state))%>%arrange(Idind,wave)%>%
 group_by(Idind)%>%mutate(treatment_weight_cum=cumprod(tp_num/tp_den))%>%ungroup()
tw<-outcome_rows%>%group_by(Idind)%>%summarise(treatment_weight=last(treatment_weight_cum),.groups="drop")

# Sequential retention weights: model observation at the next scheduled wave among those observed now.
firsts<-d%>%filter(!is.na(fuel_state))%>%group_by(Idind)%>%summarise(first_wave=min(wave),.groups="drop")
grid<-firsts%>%rowwise()%>%reframe(Idind=Idind,wave=waves[waves>=first_wave])%>%ungroup()%>%
 left_join(d%>%select(-previous_fuel),by=c("Idind","wave"))%>%arrange(Idind,wave)%>%group_by(Idind)%>%
 mutate(obs=as.integer(!is.na(fuel_state)),prior_all_observed=cumprod(obs),next_obs=lead(obs),current_state=factor(fuel_state,levels=c("solid_only","mixed","clean_only")))%>%ungroup()
risk<-grid%>%filter(wave<2009,obs==1,prior_all_observed==1,!is.na(next_obs))%>%mutate(previous_fuel="entry")
risk<-prep(risk)
for(v in c("baseline_age","baseline_urban","baseline_education","baseline_income","baseline_assets")){
 risk[[paste0(v,"_m")]]<-is.na(risk[[v]]); med<-median(risk[[v]],na.rm=TRUE);if(!is.finite(med))med<-0;risk[[paste0(v,"_i")]]<-ifelse(is.na(risk[[v]]),med,risk[[v]])
}
risk$cp_num<-risk$cp_den<-NA_real_
cnum<-next_obs~ns(baseline_age_i,3)+sex+province_f+baseline_urban_i+baseline_education_i+baseline_income_i+baseline_assets_i
cden<-next_obs~current_state+ns(age,3)+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
for(w in waves[waves<2009]){z<-risk%>%filter(wave==w);if(nrow(z)<200)next;fn<-glm(cnum,data=z,family=binomial());fd<-glm(cden,data=z,family=binomial());ii<-which(risk$wave==w);risk$cp_num[ii]<-pmax(pmin(predict(fn,newdata=risk[ii,],type="response"),.999),.001);risk$cp_den[ii]<-pmax(pmin(predict(fd,newdata=risk[ii,],type="response"),.999),.001)}
retw<-risk%>%filter(Idind%in%complete_ids,next_obs==1)%>%group_by(Idind)%>%summarise(retention_weight=prod(cp_num/cp_den,na.rm=TRUE),.groups="drop")

# Biomarker participation at 2009 among continuously observed participants with a valid 2009 state.
sel<-grid%>%filter(wave==2009,prior_all_observed==1,obs==1)%>%mutate(previous_fuel="entry")%>%prep()
for(v in c("baseline_age","baseline_urban","baseline_education","baseline_income","baseline_assets")){
 sel[[paste0(v,"_m")]]<-is.na(sel[[v]]);med<-median(sel[[v]],na.rm=TRUE);if(!is.finite(med))med<-0;sel[[paste0(v,"_i")]]<-ifelse(is.na(sel[[v]]),med,sel[[v]])
}
bnum<-biomarker2009~ns(baseline_age_i,3)+sex+province_f+baseline_urban_i+baseline_education_i+baseline_income_i+baseline_assets_i
bden<-biomarker2009~current_state+ns(age,3)+sex+province_f+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i+urban_m+education_m+income_m+assets_m+smoking_m+alcohol_m
fbn<-glm(bnum,data=sel,family=binomial());fbd<-glm(bden,data=sel,family=binomial())
sel$bp_num<-pmax(pmin(predict(fbn,type="response"),.999),.001);sel$bp_den<-pmax(pmin(predict(fbd,type="response"),.999),.001)
bw<-sel%>%filter(biomarker2009==1)%>%transmute(Idind,biomarker_weight=bp_num/bp_den)

weights<-support%>%filter(complete_sequence)%>%select(Idind,continued_any_solid,sustained_clean,reversal)%>%
 left_join(tw,by="Idind")%>%left_join(retw,by="Idind")%>%left_join(bw,by="Idind")%>%
 mutate(retention_weight=ifelse(is.na(retention_weight),1,retention_weight),biomarker_weight=ifelse(is.na(biomarker_weight),1,biomarker_weight),
  censoring_weight=retention_weight*biomarker_weight,combined_weight=treatment_weight*censoring_weight,
  strategy=case_when(continued_any_solid~"continued_any_solid",sustained_clean~"sustained_clean_only",TRUE~"other"))
if(any(!is.finite(weights$combined_weight)))stop("Non-finite combined weights")
saveRDS(weights,file.path(dd,"causal_participant_weights_internal_restricted.rds"),compress="gzip")
saveRDS(outcome_rows,file.path(dd,"causal_treatment_long_internal_restricted.rds"),compress="gzip")
write.csv(weights%>%select(-Idind),file.path(od,"causal_participant_weights_deidentified.csv"),row.names=FALSE)
cat("FINAL WEIGHTS N=",nrow(weights)," A=",sum(weights$strategy=="continued_any_solid")," B=",sum(weights$strategy=="sustained_clean_only"),"\n",sep="")
cat("IPTW_IPCW=PASS\n")
