options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages(library(dplyr))
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);dd<-file.path(project,"data","causal_landmark");od<-file.path(project,"results","causal_landmark");td<-file.path(project,"tables","causal_landmark");ld<-file.path(project,"logs","causal_landmark")
g<-read.csv(file.path(od,"landmark_causal_gate_by_method.csv"));fg<-read.csv(file.path(od,"landmark_final_weight_gate.csv"));lm<-readRDS(file.path(dd,"landmark_analysis_no_outcomes_internal_restricted.rds"))
fuel09<-lm%>%mutate(treatment2006=ifelse(A==1,"clean-only 2006","any-solid 2006"),fuel2009=ifelse(is.na(fuel2009),"missing",as.character(fuel2009)))%>%count(treatment2006,fuel2009,name="N")%>%group_by(treatment2006)%>%mutate(proportion=N/sum(N))%>%ungroup();write.csv(fuel09,file.path(td,"Table_landmark_2009_fuel_descriptive.csv"),row.names=FALSE)
sens<-bind_rows(g%>%transmute(analysis=ifelse(method=="glm","Primary GLM overlap treatment weight","CBPS overlap treatment-weight sensitivity"),balance_pass,ESS_pass,overall_pass,interpretation=ifelse(overall_pass,"Treatment-only gate passed","Treatment-only gate failed")),data.frame(analysis="GLM overlap x biomarker IPCW",balance_pass=fg$balance_pass,ESS_pass=fg$ESS_pass,overall_pass=FALSE,interpretation="Final gate failed; no outcome sensitivity models permitted"))
write.csv(sens,file.path(td,"Table_landmark_sensitivity.csv"),row.names=FALSE);writeLines("LANDMARK_SENSITIVITY=DIAGNOSTICS_ONLY",file.path(ld,"08_landmark_sensitivity.log"),useBytes=TRUE)

