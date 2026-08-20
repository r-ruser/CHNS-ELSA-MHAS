options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork);library(scales);library(ragg)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
base<-file.path(project,"causal_dml");rd<-file.path(base,"results");fd<-file.path(base,"figures");dir.create(fd,recursive=TRUE,showWarnings=FALSE)
blue<-"#0072B2";orange<-"#D55E00";grey<-"#6B6B6B"
theme_nature<-theme_classic(base_size=9,base_family="Arial")+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=8,color=grey),axis.title=element_text(size=8.5),axis.text=element_text(size=8),strip.background=element_blank(),strip.text=element_text(face="bold",size=8.5),legend.position="none",plot.margin=margin(5,7,5,7))

eff<-read.csv(file.path(rd,"Table_DML_primary_effects.csv"))%>%mutate(outcome=recode(outcome,HD_z="HD (SD)",KDM_BAA="KDM-BAA (years)"))
pa<-ggplot(eff,aes(estimate,1))+geom_vline(xintercept=0,lty=2,color=grey,linewidth=.4)+geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation="y",width=0,color=blue,linewidth=.75)+geom_point(shape=21,fill="white",color=blue,size=2.7,stroke=.8)+facet_wrap(~outcome,scales="free_x",nrow=1)+labs(title="a  Cross-fitted Super Learner DML",subtitle="Clean-only vs continued any-solid; overlap population",x="Marginal mean difference (95% CI)",y=NULL)+theme_nature+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

diag<-bind_rows(lapply(1:5,function(j)read.csv(file.path(rd,sprintf("dml_diagnostics_imp%d.csv",j)))))%>%group_by(imputation)%>%slice(1)%>%ungroup()%>%transmute(imputation,`Clean combined ESS`=ess_clean_observed_combined,`Any-solid combined ESS`=ess_anysolid_observed_combined)%>%tidyr::pivot_longer(-imputation,names_to="group",values_to="ESS")
pb<-ggplot(diag,aes(factor(imputation),ESS,color=group,group=group))+geom_line(linewidth=.7)+geom_point(size=2)+scale_color_manual(values=c(blue,orange),labels=c("Any-solid","Clean-only"))+labs(title="b  Effective sample size",subtitle="Overlap plus 2009 observation correction",x="MICE dataset",y="ESS",color=NULL)+theme_nature+theme(legend.position="bottom",legend.justification="left")

sh<-read.csv(file.path(rd,"SHAP_importance_all_outcomes.csv"))%>%group_by(outcome)%>%slice_max(mean_abs_SHAP,n=8,with_ties=FALSE)%>%ungroup()%>%mutate(outcome=recode(outcome,HD_z="HD",KDM_BAA="KDM-BAA"),feature=gsub("_"," ",feature),feature=reorder(feature,mean_abs_SHAP))
pc<-ggplot(sh,aes(mean_abs_SHAP,feature,fill=outcome))+geom_col(width=.7)+facet_wrap(~outcome,scales="free",nrow=1)+scale_fill_manual(values=c(blue,orange))+labs(title="c  SHAP importance",subtitle="Conditional-contrast prediction; not a causal decomposition",x="Mean |SHAP value|",y=NULL)+theme_nature

shorten<-function(x){x<-gsub("proportion prior","prop. prior",x);x<-gsub("province","prov.",x);x<-gsub("last two wave pattern","last-2-wave",x);x<-gsub("ever clean before 2006","prior clean",x);x<-gsub("ever mixed before 2006","prior mixed",x);x<-gsub("duration weighted solid equivalent history to 2004","solid-history",x);substr(x,1,34)}
iq<-bind_rows(lapply(c("HD_z","KDM_BAA"),function(y)read.csv(file.path(rd,paste0("SHAPIQ_interactions_",y,".csv")))))%>%mutate(f1=shorten(gsub("_"," ",feature_1)),f2=shorten(gsub("_"," ",feature_2)),pair=paste(f1,f2,sep=" x "))%>%group_by(outcome,pair)%>%summarise(mean_abs_kSII=mean(abs(kSII)),signed=mean(kSII),.groups="drop")%>%group_by(outcome)%>%slice_max(mean_abs_kSII,n=6,with_ties=FALSE)%>%ungroup()%>%mutate(outcome=recode(outcome,HD_z="HD",KDM_BAA="KDM-BAA"),pair=reorder(pair,mean_abs_kSII))
pd<-ggplot(iq,aes(mean_abs_kSII,pair,fill=signed>=0))+geom_col(width=.7)+facet_wrap(~outcome,scales="free",nrow=1)+scale_fill_manual(values=c(`TRUE`=orange,`FALSE`=blue))+scale_x_continuous(breaks=pretty_breaks(n=3))+labs(title="d  SHAP-IQ pairwise interactions",subtitle="Top-eight conditional slice; predictive interaction only (orange positive, blue negative)",x="Mean |k-SII|",y=NULL)+theme_nature

fig<-(pa|pb)/pc/pd+plot_layout(heights=c(.9,1.15,1.15))+plot_annotation(caption="ATO: average treatment effect in the overlap population. DML uses five community-grouped outer folds within each of five MICE datasets.",theme=theme(plot.caption=element_text(size=7,color=grey,hjust=0)))
stem<-file.path(fd,"Figure_DML_SHAP_SHAPIQ")
ggsave(paste0(stem,".pdf"),fig,width=183,height=240,units="mm",device=cairo_pdf)
ggsave(paste0(stem,".svg"),fig,width=183,height=240,units="mm",device=svglite::svglite)
agg_png(paste0(stem,".png"),width=2161,height=2835,units="px",res=300);print(fig);dev.off()
agg_tiff(paste0(stem,".tiff"),width=4323,height=5669,units="px",res=600,compression="lzw");print(fig);dev.off()
writeLines(c("NATURE_COMPOSITE=PASS","R-only composite; 183 x 240 mm; TIFF 600 dpi; source data are the panel CSV files."),file.path(base,"logs","24_nature_composite.log"),useBytes=TRUE)
