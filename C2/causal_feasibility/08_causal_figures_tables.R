options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(cowplot);library(grid);library(svglite);library(ragg)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");fd<-file.path(project,"figures","causal");ld<-file.path(project,"logs","causal");dir.create(fd,FALSE,TRUE)
bal<-read.csv(file.path(od,"balance_wave_summary.csv"));wd<-read.csv(file.path(td,"Supplementary_Table_weight_diagnostics.csv"))%>%filter(weight=="combined")
bal_methods<-bind_rows(
 bal%>%mutate(method="Multinomial"),
 read.csv(file.path(od,"CBPS_balance_wave_summary.csv"))%>%mutate(method="CBPS"),
 read.csv(file.path(od,"overlap_balance_wave_summary.csv"))%>%mutate(method="Overlap"),
 read.csv(file.path(od,"SL_balance_wave_summary.csv"))%>%mutate(method="SL-IPTW"))
ess_methods<-bind_rows(
 data.frame(method="Multinomial",minimum_strategy_ESS=min(read.csv(file.path(od,"primary_weight_group_ESS.csv"))$ESS)),
 data.frame(method="CBPS",minimum_strategy_ESS=min(read.csv(file.path(od,"CBPS_primary_group_ESS.csv"))$ESS)),
 data.frame(method="Overlap",minimum_strategy_ESS=min(read.csv(file.path(od,"overlap_primary_group_ESS.csv"))$ESS_treatment)),
 data.frame(method="SL-IPTW",minimum_strategy_ESS=min(read.csv(file.path(od,"SL_primary_group_ESS.csv"))$ESS)))
theme_n<-theme_classic(base_size=8.2,base_family="sans")+theme(plot.title=element_text(face="bold",size=9.5),plot.subtitle=element_text(size=7.2,colour="#666666"),axis.title=element_text(size=8),axis.text=element_text(size=7),legend.position="none",plot.margin=margin(5,5,5,5))

# Panel A: vector DAG drawn natively in R.
nodes<-data.frame(name=c("Baseline\nconfounders","SES(t)","Fuel(t)","SES(t+1)","Fuel(t+1)","Retention /\nbiomarker","HD 2009"),x=c(0,1.2,1.2,2.5,2.5,3.8,5),y=c(1.5,2.25,.75,2.25,.75,1.5,1.5),type=c("base","ses","fuel","ses","fuel","selection","outcome"))
edges<-data.frame(from=c(1,1,1,2,2,3,3,3,4,4,5,2,6),to=c(2,3,7,3,4,4,5,6,5,7,7,6,7))%>%left_join(nodes%>%mutate(from=row_number())%>%select(from,x,y),by="from")%>%left_join(nodes%>%mutate(to=row_number())%>%select(to,x2=x,y2=y),by="to")
cols<-c(base="#BDBDBD",ses="#009E73",fuel="#0072B2",selection="#8E6BBE",outcome="#D55E00")
pA<-ggplot()+geom_segment(data=edges,aes(x=x+.18,y=y,xend=x2-.2,yend=y2),arrow=arrow(length=unit(1.5,"mm"),type="closed"),linewidth=.35,colour="#555555")+geom_label(data=nodes,aes(x,y,label=name,fill=type),size=2.4,linewidth=.25,label.r=unit(1.2,"mm"),colour="white",fontface="bold")+scale_fill_manual(values=cols)+coord_cartesian(xlim=c(-.4,5.45),ylim=c(.2,2.75),clip="off")+labs(title="a  Longitudinal causal structure",subtitle="Fuel(t) to SES(t+1) to Fuel(t+1) creates treatment-confounder feedback")+theme_void(base_family="sans")+theme(plot.title=element_text(face="bold",size=9.5),plot.subtitle=element_text(size=7.2,colour="#666666"),legend.position="none",plot.margin=margin(8,5,3,5))

pB<-ggplot(bal_methods%>%filter(is.finite(max_abs_SMD)),aes(wave,max_abs_SMD,colour=method,group=method))+geom_hline(yintercept=.1,linetype=2,colour="#B2182B",linewidth=.5)+geom_line(linewidth=.65)+geom_point(size=1.45)+scale_colour_manual(values=c("Multinomial"="#E69F00","CBPS"="#0072B2","Overlap"="#CC79A7","SL-IPTW"="#009E73"))+scale_x_continuous(breaks=bal$wave)+scale_y_continuous(expand=expansion(mult=c(0,.03)))+labs(title="b  No weighting approach achieved balance",subtitle="Maximum pairwise absolute SMD; dashed line marks 0.10",x="Treatment wave",y="Maximum |SMD|",colour=NULL)+guides(colour=guide_legend(nrow=1))+theme_n+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="top",legend.text=element_text(size=6.3),legend.key.width=unit(5,"mm"))
pC<-ggplot(ess_methods,aes(reorder(method,minimum_strategy_ESS),minimum_strategy_ESS,fill=method))+geom_hline(yintercept=50,linetype=2,colour="#B2182B",linewidth=.5)+geom_col(width=.62)+geom_text(aes(label=round(minimum_strategy_ESS)),hjust=-.2,size=2.5)+coord_flip()+scale_fill_manual(values=c("Multinomial"="#E69F00","CBPS"="#0072B2","Overlap"="#CC79A7","SL-IPTW"="#009E73"))+expand_limits(y=max(ess_methods$minimum_strategy_ESS)*1.14)+labs(title="c  Overlap support collapsed",subtitle="Minimum strategy ESS; dashed line marks 50",x=NULL,y="Effective sample size")+theme_n+theme(axis.text.y=element_text(size=6.6))
pD<-ggplot()+annotate("rect",xmin=0,xmax=1,ymin=0,ymax=1,fill="#F7F7F7",colour="#B2182B",linewidth=.7)+annotate("text",x=.025,y=.78,hjust=0,label="d  Prespecified causal stop rule",fontface="bold",size=3.2,family="sans")+annotate("text",x=.025,y=.48,hjust=0,label="Multinomial IPTW, CBPS, generalized overlap and Super Learner-IPTW all failed longitudinal balance.",size=2.7,family="sans")+annotate("text",x=.025,y=.20,hjust=0,label="MSM effects: WITHHELD     |     g-formula interventions: WITHHELD",fontface="bold",colour="#B2182B",size=2.9,family="sans")+annotate("text",x=.975,y=.78,hjust=1,label="Descriptive multistate and GAM results remain unchanged.",size=2.5,family="sans",colour="#555555")+coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off")+theme_void()
fig<-plot_grid(pA,plot_grid(pB,pC,ncol=2,rel_widths=c(1,1)),pD,ncol=1,rel_heights=c(.72,1.15,.42))
svglite::svglite(file.path(fd,"Figure3_longitudinal_causal_framework.svg"),width=183/25.4,height=190/25.4,system_fonts=list(sans="Arial"));print(fig);dev.off()
cairo_pdf(file.path(fd,"Figure3_longitudinal_causal_framework.pdf"),width=183/25.4,height=190/25.4,family="sans");print(fig);dev.off()
ragg::agg_tiff(file.path(fd,"Figure3_longitudinal_causal_framework.tiff"),width=4323,height=4488,units="px",res=600,compression="lzw");print(fig);dev.off()
ragg::agg_png(file.path(fd,"Figure3_longitudinal_causal_framework.png"),width=2161,height=2244,units="px",res=300);print(fig);dev.off()
write.csv(bal_methods,file.path(fd,"Figure3_source_data_panel_b.csv"),row.names=FALSE);write.csv(ess_methods,file.path(fd,"Figure3_source_data_panel_c.csv"),row.names=FALSE)
writeLines(c("Figure 3 | Longitudinal causal framework and identification diagnostics.","a, Longitudinal DAG. b, Maximum absolute standardized difference after multinomial IPTW, CBPS-IPTW, generalized overlap weighting and Super Learner-IPTW at each treatment wave; dashed line marks 0.10. c, Minimum strategy-specific ESS for each approach; generalized overlap uses treatment-only ESS because longitudinal overlap support is the estimand-defining quantity. d, The prespecified balance stop rule was triggered, so MSM and g-formula effects were not estimated. SMD, standardized mean difference; ESS, effective sample size; IPTW, inverse probability of treatment weighting."),file.path(fd,"Figure3_legend.md"),useBytes=TRUE)
writeLines(c("R-only export: PASS","Canvas: 183 x 190 mm","TIFF: 4323 x 4488 px at 600 dpi","PNG: 2161 x 2244 px at 300 dpi","Vector SVG/PDF: PASS","Visual QA pending scripted preview inspection."),file.path(fd,"Figure3_qa_notes.md"),useBytes=TRUE)
writeLines("CAUSAL_FIGURES_TABLES=PASS",file.path(ld,"08_causal_figures_tables.log"),useBytes=TRUE)
