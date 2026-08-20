options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages(library(dplyr))
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal")
wd<-read.csv(file.path(td,"Supplementary_Table_weight_diagnostics.csv"));bal<-read.csv(file.path(od,"balance_wave_summary.csv"));ess<-read.csv(file.path(od,"primary_weight_group_ESS.csv"))
sens<-bind_rows(
 wd%>%filter(weight=="combined")%>%transmute(analysis=paste0("Combined weights: ",rule),diagnostic="ESS",value=ESS,threshold="overall ESS reported",status=ifelse(ESS>=100,"stable overall","unstable")),
 bal%>%transmute(analysis=paste0("Wave ",wave),diagnostic="maximum absolute SMD",value=max_abs_SMD,threshold="<0.10",status=ifelse(!is.na(value)&value<.10,"pass","fail")),
 ess%>%transmute(analysis=paste0("Primary 1/99: ",strategy),diagnostic="group ESS",value=ESS,threshold=">=50",status=ifelse(ESS>=50,"pass","fail")))
write.csv(sens,file.path(td,"Supplementary_Table_causal_sensitivity.csv"),row.names=FALSE)
writeLines(c("CAUSAL_SENSITIVITY=COMPLETED","Untruncated and 1/99, 2.5/97.5, 5/95 weight tails were audited.","A nonlinear prior-treatment interaction specification did not resolve wave-specific imbalance.","Outcome-based sensitivity analyses were not run after the identification stop rule."),file.path(ld,"07_causal_sensitivity.log"),useBytes=TRUE)

