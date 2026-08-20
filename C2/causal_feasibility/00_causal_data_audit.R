options(stringsAsFactors=FALSE,warn=1,cli.unicode=FALSE)
suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr);library(readr)})
set.seed(20260815)
project<-normalizePath(file.path(getwd()),winslash="/",mustWork=TRUE)
if(basename(project)!="longitudinal_upgrade")stop("Run from C2/longitudinal_upgrade")
c2<-dirname(project)
od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal")
ld<-file.path(project,"logs","causal");dd<-file.path(project,"data","causal")
for(z in c(od,td,ld,dd,file.path(project,"figures","causal")))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"00_causal_data_audit.log"),"wt",encoding="UTF-8")
sink(logcon,split=TRUE);sink(logcon,type="message")
on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
cat("Causal data audit",format(Sys.time()),"\n",R.version.string,"\n")

raw_root<-Sys.getenv("CHNS_RAW_ROOT", unset="")
if(!nzchar(raw_root) || !dir.exists(raw_root))stop("Set CHNS_RAW_ROOT to the local CHNS-STATA directory before running this script")
file_map<-c(asset="asset_12.dta",surveys="surveys_pub_12.dta",master="mast_pub_12.dta",
 pexam="pexam_pub_12.dta",education="educ_12.dta",hh_income="hhinc_10.dta")
path_for<-function(x){z<-list.files(raw_root,pattern=paste0("^",file_map[[x]],"$"),recursive=TRUE,full.names=TRUE,ignore.case=TRUE);if(length(z)!=1)stop("Expected one mapped input for ",x,"; found ",length(z));z}
waves<-c(1989,1991,1993,1997,2000,2004,2006,2009)

asset<-read_dta(path_for("asset"))%>%filter(wave%in%waves)
survey<-read_dta(path_for("surveys"))%>%filter(wave%in%waves)%>%
 select(Idind,hhid,wave,age,commid,province=t1,stratum,urban)
master<-read_dta(path_for("master"))%>%select(Idind,gender)
pexam<-read_dta(path_for("pexam"))%>%filter(wave%in%waves)%>%
 select(IDind,wave,any_of(c("U25","U40")))%>%rename(Idind=IDind)
educ<-read_dta(path_for("education"))%>%filter(wave%in%waves)%>%
 select(IDind,wave,any_of(c("A11","A12")))%>%rename(Idind=IDind)
income<-read_dta(path_for("hh_income"))%>%filter(wave%in%waves)%>%select(hhid,wave,hhincpc_cpi)

asset_vars<-intersect(c("L23","L27","L31","L105","L110","L115","L120","L140E"),names(asset))
asset2<-asset%>%mutate(across(all_of(asset_vars),~ifelse(.x%in%c(0,1),.x,NA_real_)))
asset2$asset_count<-apply(as.data.frame(asset2[,asset_vars,drop=FALSE]),1,function(z)if(all(is.na(z)))NA_real_ else sum(z,na.rm=TRUE))
asset2<-asset2%>%select(hhid,wave,L8_1,L8_2,asset_count)

fuel_class<-function(x)case_when(x%in%c(1,6,7)~"solid",x%in%c(2,4,5)~"clean",TRUE~NA_character_)
outcome<-readRDS(file.path(project,"data","dynamic_analysis_internal_restricted.rds"))%>%
 select(Idind,HD_z,KDM_BAA)%>%distinct()
panel<-survey%>%left_join(master,by="Idind")%>%left_join(asset2,by=c("hhid","wave"))%>%
 left_join(pexam,by=c("Idind","wave"))%>%left_join(educ,by=c("Idind","wave"))%>%
 left_join(income,by=c("hhid","wave"))%>%
 mutate(primary=fuel_class(L8_1),secondary=fuel_class(L8_2),
  fuel_state=case_when((primary=="solid"|secondary=="solid")&(primary=="clean"|secondary=="clean")~"mixed",
   (primary=="solid"|secondary=="solid")~"solid_only",
   (primary=="clean"|secondary=="clean")~"clean_only",TRUE~NA_character_),
  fuel_state=factor(fuel_state,levels=c("solid_only","mixed","clean_only")),
  sex=factor(gender,levels=c(1,2),labels=c("Male","Female")),
  smoking=case_when(U25%in%c(0,1)~as.numeric(U25),TRUE~NA_real_),
  alcohol=case_when(U40%in%c(0,1)~as.numeric(U40),TRUE~NA_real_),
  education_years=case_when(A11>=0~as.numeric(A11),TRUE~NA_real_),
  education_level=case_when(A12%in%0:6~as.numeric(A12),TRUE~NA_real_),
  income_asinh=asinh(hhincpc_cpi),biomarker2009=as.integer(Idind%in%outcome$Idind))%>%
 arrange(Idind,wave)%>%group_by(Idind)%>%mutate(previous_fuel=lag(fuel_state),
  next_scheduled=lead(wave),retained_next=as.integer(!is.na(lead(wave))))%>%ungroup()%>%
 left_join(outcome,by="Idind")
stopifnot(!anyDuplicated(panel[c("Idind","wave")]))

# Causal risk set: adult at first valid fuel observation, at least two valid fuel waves.
eligible_ids<-panel%>%filter(age>=18,!is.na(fuel_state))%>%group_by(Idind)%>%
 summarise(n_valid=n(),first_wave=min(wave),.groups="drop")%>%filter(n_valid>=2)%>%pull(Idind)
long<-panel%>%filter(Idind%in%eligible_ids)%>%group_by(Idind)%>%arrange(wave,.by_group=TRUE)%>%
 mutate(baseline_age=first(age[!is.na(age)]),baseline_urban=first(urban[!is.na(urban)]),
  baseline_education=first(education_years[!is.na(education_years)]),
  baseline_income=first(income_asinh[!is.na(income_asinh)]),baseline_assets=first(asset_count[!is.na(asset_count)]))%>%ungroup()

wave_audit<-long%>%group_by(wave)%>%summarise(n=n(),ids=n_distinct(Idind),fuel_known=sum(!is.na(fuel_state)),
 biomarker_ids=sum(biomarker2009==1),.groups="drop")
write.csv(wave_audit,file.path(od,"causal_wave_audit.csv"),row.names=FALSE)
saveRDS(long,file.path(dd,"causal_long_full_internal_restricted.rds"),compress="gzip")

audit_vars<-data.frame(
 variable=c("participant_id","community_id","survey_year","fuel_state","previous_fuel_state","age","sex","province","urbanicity","education","income","household_assets","smoking","alcohol","retention","biomarker2009","HD","KDM_BAA"),
 actual_CHNS_name=c("Idind","commid","wave","L8_1+L8_2","lag(L8_1+L8_2)","age","gender","t1","urban","A11/A12","hhincpc_cpi",paste(asset_vars,collapse="+"),"U25","U40","lead(observed wave)","biomarker_09 linkage","HD_z","KDM_BAA"),
 role=c("identifier","clustering","time","time-varying exposure","treatment history","baseline/time-varying confounder","baseline fixed confounder","baseline/time-varying confounder","time-varying confounder","time-varying confounder/possible mediator","time-varying confounder/possible mediator","time-varying confounder/possible mediator","time-varying behavioural confounder","time-varying behavioural confounder","censoring predictor/outcome","selection outcome","outcome component","outcome component"),
 stringsAsFactors=FALSE)
lookup<-list(age=long$age,sex=long$sex,province=long$province,urbanicity=long$urban,education=long$education_years,
 income=long$income_asinh,household_assets=long$asset_count,smoking=long$smoking,alcohol=long$alcohol,
 fuel_state=long$fuel_state,previous_fuel_state=long$previous_fuel,retention=long$retained_next,biomarker2009=long$biomarker2009,HD=long$HD_z,KDM_BAA=long$KDM_BAA)
audit_vars$waves_available<-vapply(audit_vars$variable,function(v){if(v%in%names(lookup))paste(sort(unique(long$wave[!is.na(lookup[[v]])])),collapse=";")else paste(waves,collapse=";")},character(1))
audit_vars$missingness<-vapply(audit_vars$variable,function(v){if(v%in%names(lookup))mean(is.na(lookup[[v]])) else 0},numeric(1))
audit_vars$time_varying<-audit_vars$variable%in%c("fuel_state","previous_fuel_state","age","province","urbanicity","education","income","household_assets","smoking","alcohol","retention")
audit_vars$included_in_DAG<-!audit_vars$variable%in%c("participant_id","community_id","survey_year","previous_fuel_state","KDM_BAA")
audit_vars$included_in_treatment_model<-audit_vars$variable%in%c("fuel_state","previous_fuel_state","age","sex","province","urbanicity","education","income","household_assets","smoking","alcohol")
audit_vars$included_in_censoring_model<-audit_vars$variable%in%c("fuel_state","age","sex","province","urbanicity","education","income","household_assets","smoking","alcohol","retention","biomarker2009")
audit_vars$reason<-c("linkage only","cluster-robust inference","temporal index","multinomial treatment","treatment history","confounding","baseline confounding","geography/confounding","time-varying SES","time-varying SES and possible pathway","time-varying SES and possible pathway","time-varying SES and possible pathway","behavioural confounding","behavioural confounding","loss to follow-up","2009 selection","primary outcome","secondary outcome")
write.csv(audit_vars,file.path(td,"causal_variable_audit.csv"),row.names=FALSE,fileEncoding="UTF-8")
cat("WAVES:",paste(sort(unique(long$wave)),collapse=","),"\n")
cat("ELIGIBLE IDS:",n_distinct(long$Idind)," LONG ROWS:",nrow(long)," BIOMARKER IDS:",n_distinct(long$Idind[long$biomarker2009==1]),"\n")
cat("CAUSAL_DATA_AUDIT=PASS\n")
