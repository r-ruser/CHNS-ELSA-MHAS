options(stringsAsFactors=FALSE,warn=1)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal")
status<-data.frame(method="parametric g-formula",natural_course_validation="NOT ATTEMPTED",intervention_estimate="WITHHELD",reason="The common longitudinal identification set failed the treatment-balance stop rule; policy simulation would not repair lack of exchangeability/overlap.")
write.csv(status,file.path(td,"Supplementary_Table_gformula_interventions.csv"),row.names=FALSE)
write.csv(status,file.path(od,"gformula_status.csv"),row.names=FALSE)
writeLines(c("GFORMULA=NOT_RUN","Reason: upstream identification hard-stop; no intervention simulation was interpreted."),file.path(ld,"06_gformula.log"),useBytes=TRUE)

