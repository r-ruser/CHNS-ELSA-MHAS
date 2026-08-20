options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");ld<-file.path(project,"logs","causal")
for(z in c(od,td,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"02_target_trial_specification.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);on.exit({sink();close(logcon)},add=TRUE)
d<-readRDS(file.path(dd,"causal_long_full_internal_restricted.rds"))%>%filter(biomarker2009==1,!is.na(fuel_state))%>%arrange(Idind,wave)
waves<-c(1989,1991,1993,1997,2000,2004,2006,2009)
support<-d%>%group_by(Idind)%>%summarise(first_wave=min(wave),last_wave=max(wave),n_waves=n(),person_waves=n(),
 complete_sequence=all(waves[waves>=first_wave]%in%wave),
 continued_any_solid=all(fuel_state!="clean_only"),
 sustained_clean={x<-as.character(fuel_state);k<-which(x=="clean_only");length(k)>0&&any(x[seq_len(min(k))]!="clean_only")&&all(x[min(k):length(x)]=="clean_only")},
 reversal={x<-as.character(fuel_state);k<-which(x=="clean_only");length(k)>0&&any(seq_along(x)>min(k)&x!="clean_only")},.groups="drop")
support_counts<-data.frame(strategy=c("A continued any-solid","B transition then sustained clean-only","C clean-only followed by reversal"),
 participants=c(sum(support$continued_any_solid),sum(support$sustained_clean),sum(support$reversal)),
 person_waves=c(sum(support$person_waves[support$continued_any_solid]),sum(support$person_waves[support$sustained_clean]),sum(support$person_waves[support$reversal])))
write.csv(support_counts,file.path(od,"target_strategy_support.csv"),row.names=FALSE)
complete_counts<-support%>%filter(complete_sequence)%>%summarise(complete_n=n(),A=sum(continued_any_solid),B=sum(sustained_clean),C=sum(reversal))
write.csv(complete_counts,file.path(od,"target_strategy_support_complete_sequence.csv"),row.names=FALSE)
final_contrast<-if(complete_counts$A>=200&&complete_counts$B>=200)"Strategy B versus Strategy A" else "NOT IDENTIFIABLE: insufficient sustained-strategy support"
spec<-data.frame(element=c("Eligibility criteria","Treatment strategies","Assignment procedure","Start of follow-up","End of follow-up","Outcome","Causal contrast","Estimand","Analysis method"),
 specification=c("Adult at first valid fuel observation; at least two valid fuel waves; valid 2009 biomarker outcome for the outcome analysis",
 "A: continued any-solid use (solid-only or mixed at every observed wave); B: transition from any-solid to clean-only followed by sustained clean-only use",
 "Observational emulation using sequential stabilized treatment and censoring weights","First eligible observed CHNS fuel wave","2009 biomarker examination",
 "Primary: HD_z; secondary: cross-fitted KDM_BAA in years",final_contrast,
 "E[HD2009^B]-E[HD2009^A] and E[KDM_BAA2009^B]-E[KDM_BAA2009^A] in the eligible linked cohort",
 "Marginal structural mean model with stabilized multinomial IPTW and IPCW; community-clustered robust inference"))
write.csv(spec,file.path(td,"Table1_target_trial_specification.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(c("FINAL CAUSAL ESTIMAND",spec$specification[spec$element=="Estimand"],"","POSITIVITY SUPPORT",capture.output(print(support_counts)),"",paste("Complete-sequence support:",paste(names(complete_counts),complete_counts,collapse="; "))),file.path(od,"target_trial_estimand.txt"))
saveRDS(support,file.path(dd,"target_strategy_support_internal_restricted.rds"),compress="gzip")
cat("FINAL CONTRAST:",final_contrast,"\n");print(support_counts);print(complete_counts);cat("TARGET_TRIAL_SPEC=PASS\n")
