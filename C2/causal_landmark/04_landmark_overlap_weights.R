options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(WeightIt);library(cobalt);library(ggplot2);library(svglite)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);dd<-file.path(project,"data","causal_landmark");od<-file.path(project,"results","causal_landmark");td<-file.path(project,"tables","causal_landmark");fd<-file.path(project,"figures","causal_landmark");ld<-file.path(project,"logs","causal_landmark")
logcon<-file(file.path(ld,"04_landmark_overlap_weights.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
ims<-readRDS(file.path(dd,"landmark_MICE5_no_outcomes_internal_restricted.rds"))
form<-A~age2004+sex+province+urban_recent+education_recent+income_recent+income_mean+income_change+assets_recent+assets_mean+assets_change+smoking_recent+alcohol_recent+fuel2004+number_prior_fuel_waves+proportion_prior_solid_only+proportion_prior_mixed+ever_clean_before_2006+ever_mixed_before_2006+number_fuel_transitions+number_clean_to_mixed_reversals+number_mixed_to_solid_reversals+duration_weighted_solid_equivalent_history_to_2004+last_two_wave_pattern
weights<-list();balance<-list();psdiag<-list();ess<-list();cbps_status<-list()
for(j in 1:5){
 z<-ims[[j]];z$A<-as.integer(z$A);z$A_f<-factor(z$A,levels=c(0,1),labels=c("continued_any_solid","clean_only"))
 ow<-WeightIt::weightit(form,data=z,method="glm",estimand="ATO");z$ps_glm<-ow$ps;z$w_glm<-ow$weights
 cb<-tryCatch(WeightIt::weightit(form,data=z,method="cbps",estimand="ATO",over=TRUE),error=function(e)tryCatch(WeightIt::weightit(form,data=z,method="cbps",estimand="ATO",over=FALSE),error=function(e2)NULL))
 if(!is.null(cb)){z$ps_cbps<-cb$ps;z$w_cbps<-cb$weights;cbps_status[[j]]<-data.frame(imputation=j,estimable=TRUE)}else{z$ps_cbps<-z$w_cbps<-NA_real_;cbps_status[[j]]<-data.frame(imputation=j,estimable=FALSE)}
 for(method in c("glm","cbps")){w<-z[[paste0("w_",method)]];if(all(is.na(w)))next
  bt<-cobalt::bal.tab(form,data=z,weights=w,method="weighting",estimand="ATO",un=TRUE,quick=FALSE)
  b<-as.data.frame(bt$Balance);b$covariate<-rownames(b);rownames(b)<-NULL;b$imputation<-j;b$method<-method;balance[[length(balance)+1]]<-b
  e<-z[[paste0("ps_",method)]];pe<-data.frame(A_f=z$A_f,e=e);psdiag[[length(psdiag)+1]]<-pe%>%group_by(A_f)%>%summarise(N=n(),minimum=min(e),P1=quantile(e,.01),P5=quantile(e,.05),P25=quantile(e,.25),median=median(e),P75=quantile(e,.75),P95=quantile(e,.95),P99=quantile(e,.99),maximum=max(e),frac_lt_001=mean(e<.01),frac_lt_005=mean(e<.05),frac_gt_095=mean(e>.95),frac_gt_099=mean(e>.99),.groups="drop")%>%mutate(imputation=j,method=method)
  wz<-data.frame(A_f=z$A_f,w=w)
  ess[[length(ess)+1]]<-bind_rows(wz%>%group_by(A_f)%>%summarise(N=n(),sum_w=sum(w),ESS=sum(w)^2/sum(w^2),.groups="drop")%>%mutate(scope="treatment-specific"),wz%>%summarise(N=n(),sum_w=sum(w),ESS=sum(w)^2/sum(w^2))%>%mutate(A_f="overall",scope="overall"))%>%mutate(imputation=j,method=method)
 }
 weights[[j]]<-z;cat("LANDMARK_WEIGHT_DONE=",j,"/5\n",sep="")
}
saveRDS(weights,file.path(dd,"landmark_overlap_weights_MICE5_internal_restricted.rds"),compress="gzip")
bal<-bind_rows(balance);psd<-bind_rows(psdiag);es<-bind_rows(ess);write.csv(bal,file.path(td,"Table_landmark_balance.csv"),row.names=FALSE);write.csv(psd,file.path(od,"landmark_propensity_diagnostics.csv"),row.names=FALSE);write.csv(es,file.path(od,"landmark_weight_ESS.csv"),row.names=FALSE);write.csv(bind_rows(cbps_status),file.path(od,"landmark_CBPS_estimability.csv"),row.names=FALSE)
# Main diagnostic plots use imputation 1; all-imputation numerical gates are separate.
z<-weights[[1]]
p1<-ggplot(z,aes(ps_glm,fill=A_f,colour=A_f))+geom_density(alpha=.18,linewidth=.65)+scale_fill_manual(values=c("continued_any_solid"="#D55E00","clean_only"="#0072B2"),labels=c("Continued any-solid","Clean-only"))+scale_colour_manual(values=c("continued_any_solid"="#D55E00","clean_only"="#0072B2"),labels=c("Continued any-solid","Clean-only"))+theme_classic(base_size=9)+labs(title="2006 treatment propensity overlap",subtitle="Imputation 1; GLM overlap-weighting propensity score",x="Pr(clean-only in 2006 | pretreatment history)",y="Density",fill=NULL,colour=NULL)+theme(legend.position="top")
ggsave(file.path(fd,"Figure_landmark_propensity_overlap.pdf"),p1,width=170,height=105,units="mm",device=cairo_pdf);ggsave(file.path(fd,"Figure_landmark_propensity_overlap.png"),p1,width=2008,height=1240,units="px",res=300,device=ragg::agg_png);ggsave(file.path(fd,"Figure_landmark_propensity_overlap.svg"),p1,width=170,height=105,units="mm",device=svglite::svglite)
bt1<-cobalt::bal.tab(form,data=z,weights=z$w_glm,method="weighting",estimand="ATO",un=TRUE,quick=FALSE);lp<-cobalt::love.plot(bt1,abs=TRUE,thresholds=c(m=.1),var.order="unadjusted",colors=c("#777777","#0072B2"),shapes=c(16,17),sample.names=c("Unweighted","Overlap weighted"),stars="raw")+theme_minimal(base_size=8)+theme(legend.position="top",panel.grid.minor=element_blank())
ggsave(file.path(fd,"Figure_landmark_loveplot.pdf"),lp,width=175,height=150,units="mm",device=cairo_pdf);ggsave(file.path(fd,"Figure_landmark_loveplot.svg"),lp,width=175,height=150,units="mm",device=svglite::svglite);ggsave(file.path(fd,"Figure_landmark_loveplot.png"),lp,width=2067,height=1772,units="px",res=300,device=ragg::agg_png)
writeLines("LANDMARK_OVERLAP_WEIGHTS=PASS_DIAGNOSTICS_ONLY",file.path(ld,"04_landmark_overlap_weights.log.status"),useBytes=TRUE)
