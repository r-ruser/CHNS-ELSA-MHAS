options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(svglite);library(ragg)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);stopifnot(basename(project)=="longitudinal_upgrade")
resdir<-file.path(project,"results");figdir<-file.path(project,"figures");dir.create(figdir,FALSE,TRUE)
q<-read.csv(file.path(resdir,"transition_intensities.csv"));occ<-read.csv(file.path(resdir,"occupancy_from_solid.csv"))
ar<-read.csv(file.path(resdir,"dynamic_association_results.csv"));bs<-read.csv(file.path(resdir,"community_bootstrap_summary.csv"));gc<-read.csv(file.path(resdir,"gam_adjusted_curves.csv"));gd<-read.csv(file.path(resdir,"gam_diagnostics.csv"))
pal<-c(solid_only="#D55E00",mixed="#E69F00",clean_only="#0072B2")
labstate<-c(solid_only="Solid only",mixed="Mixed fuels",clean_only="Clean only")
theme_nat<-theme_classic(base_size=6.5,base_family="Arial")+theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),plot.title=element_text(size=7.5,face="bold"),plot.subtitle=element_text(size=6.2,color="#555555"),axis.title=element_text(size=6.5),axis.text=element_text(size=6),strip.background=element_blank(),strip.text=element_text(size=6.5,face="bold"),legend.text=element_text(size=5.8),legend.title=element_blank(),plot.margin=margin(4,5,4,4))

# a: triangular reversible-state loop; opposite directions use separate curves.
nodes<-data.frame(state=c("solid_only","mixed","clean_only"),x=c(-.92,.92,0),y=c(-.48,-.48,.62),label=labstate[c("solid_only","mixed","clean_only")])
edges<-q%>%left_join(nodes%>%select(from=state,x,y),by="from")%>%
 left_join(nodes%>%select(to=state,xend=x,yend=y),by="to")%>%
 mutate(pair=paste(pmin(from,to),pmax(from,to),sep="__"),forward=from<to,
 label=sprintf("%s→%s  %.3f (%.3f–%.3f)",
   recode(from,solid_only="Solid",mixed="Mixed",clean_only="Clean"),
   recode(to,solid_only="solid",mixed="mixed",clean_only="clean"),
   estimate,lower,upper),
 label_x=case_when(pair=="mixed__solid_only"~0,pair=="clean_only__solid_only"~-.54,TRUE~.54),
 label_y=case_when(pair=="mixed__solid_only"&forward~-.72,pair=="mixed__solid_only"~-.88,
 pair=="clean_only__solid_only"&forward~.20,pair=="clean_only__solid_only"~.03,
 pair=="clean_only__mixed"&forward~.20,TRUE~.03),
 hjust=case_when(pair=="clean_only__solid_only"~1,pair=="clean_only__mixed"~0,TRUE~.5))
pa<-ggplot()+
 geom_curve(data=filter(edges,forward),aes(x=x,y=y,xend=xend,yend=yend,linewidth=estimate,color=from),curvature=.13,arrow=arrow(length=unit(1.7,"mm")),lineend="round",show.legend=FALSE)+
 geom_curve(data=filter(edges,!forward),aes(x=x,y=y,xend=xend,yend=yend,linewidth=estimate,color=from),curvature=.13,arrow=arrow(length=unit(1.7,"mm")),lineend="round",show.legend=FALSE)+
 geom_label(data=nodes,aes(x,y,label=label,fill=state),color="white",fontface="bold",size=2.5,linewidth=0,label.padding=unit(2.2,"mm"),show.legend=FALSE)+
 geom_text(data=edges,aes(label_x,label_y,label=label,hjust=hjust),size=1.72,color="#333333")+
 scale_color_manual(values=pal)+scale_fill_manual(values=pal)+scale_linewidth(range=c(.25,1.6))+
 coord_cartesian(xlim=c(-1.50,1.50),ylim=c(-1.02,.92),clip="off")+theme_void(base_family="Arial")+
 labs(title="a  Reversible annual fuel-state transitions",subtitle="Estimate (95% CI); arrows show direction and width scales with intensity")+theme(plot.title=element_text(size=7.5,face="bold"),plot.subtitle=element_text(size=6.2,color="#555555"),plot.margin=margin(4,8,4,4))

pb_pal<-setNames(unname(pal),unname(labstate))
pb<-occ%>%mutate(state=factor(state,levels=names(labstate),labels=labstate))%>%ggplot(aes(year,probability,color=state))+
 geom_line(linewidth=.75)+scale_color_manual(values=pb_pal)+scale_y_continuous(labels=scales::percent_format(accuracy=1),limits=c(0,1))+
 labs(title="b  Model-derived state occupancy",subtitle="Starting from solid-only fuel use",x="Years since starting state",y="Probability",color=NULL)+theme_nat+theme(legend.position=c(.74,.77),legend.background=element_rect(fill=scales::alpha("white",.85),color=NA))

term_map<-c("Solid-equivalent history"="Solid-equivalent history\n(per 10 percentage points)","Sustained clean-only duration"="Sustained clean-only duration\n(per 5 years)","Clean-to-solid/mixed reversal"="Any observed reversal\n(yes versus no)")
term_code<-c("Solid-equivalent history"="solid_prop10","Sustained clean-only duration"="clean5","Clean-to-solid/mixed reversal"="rev_any")
pcd<-ar%>%mutate(term=unname(term_code[exposure]))%>%left_join(bs,by=c("outcome","term"),suffix=c("_robust","_boot"))%>%
 mutate(label=unname(term_map[exposure]),outcome_label=ifelse(outcome=="HD_z","HD (SD)","KDM-BAA (years)"))
pc_long<-bind_rows(
 pcd%>%transmute(outcome_label,label,estimate,lower=lower_robust,upper=upper_robust,method="Community-robust",offset=.11),
 pcd%>%transmute(outcome_label,label,estimate,lower=lower_boot,upper=upper_boot,method="Cluster bootstrap",offset=-.11))%>%
 mutate(label=factor(label,levels=rev(unname(term_map))),y=as.numeric(label)+offset)
pc<-ggplot(pc_long,aes(estimate,y,color=method))+geom_vline(xintercept=0,linetype=2,color="#777777",linewidth=.35)+
 geom_errorbar(aes(xmin=lower,xmax=upper),width=.08,linewidth=.68)+geom_point(size=1.7)+
 facet_wrap(~outcome_label,scales="free_x",ncol=2)+
 scale_y_continuous(breaks=seq_along(levels(pc_long$label)),labels=levels(pc_long$label),expand=expansion(mult=c(.14,.14)))+
 scale_color_manual(values=c("Community-robust"="#0072B2","Cluster bootstrap"="#E69F00"))+
 labs(title="c  Dynamic fuel history and 2009 biological age",subtitle="Intervals are vertically separated to show agreement between inference methods",x="Adjusted mean difference",y=NULL,color=NULL)+theme_nat+theme(axis.text.y=element_text(size=5.6),legend.position="top",legend.justification="left")

gcd<-gc%>%left_join(gd[,c("outcome","edf","p_value")],by="outcome")%>%mutate(outcome_label=ifelse(outcome=="HD_z",sprintf("HD (EDF %.2f; P=%.3f)",edf,p_value),sprintf("KDM-BAA (EDF %.2f; P=%.3f)",edf,p_value)))
pd<-ggplot(gcd,aes(solid_equiv_prop,estimate))+geom_ribbon(aes(ymin=lower,ymax=upper),fill="#56B4E9",alpha=.22)+geom_line(color="#0072B2",linewidth=.7)+geom_hline(yintercept=0,linetype=2,color="#777777",linewidth=.3)+facet_wrap(~outcome_label,scales="free_y",ncol=2)+scale_x_continuous(labels=scales::percent_format(accuracy=1))+
 labs(title="d  Adjusted nonlinear sensitivity",subtitle="Cross-sectional 2009 contrast relative to 0% solid-equivalent history",x="Solid-equivalent proportion of observed history",y="Adjusted difference")+theme_nat

top<-pa+pb+plot_layout(ncol=2,widths=c(1,1))
bottom<-pc+pd+plot_layout(ncol=1,heights=c(1.08,1))
fig<-top/bottom+plot_layout(heights=c(1.02,1.42))
base<-file.path(figdir,"Figure4_C2_dynamic_fuel_biological_age");w<-183/25.4;h<-150/25.4
svglite(base%>%paste0(".svg"),width=w,height=h);print(fig);dev.off()
cairo_pdf(paste0(base,".pdf"),width=w,height=h,family="Arial");print(fig);dev.off()
agg_tiff(paste0(base,".tiff"),width=4323,height=3543,units="px",res=600,compression="lzw");print(fig);dev.off()
agg_png(paste0(base,".png"),width=2161,height=1772,units="px",res=300);print(fig);dev.off()
write.csv(edges,file.path(figdir,"Figure4_source_panel_a.csv"),row.names=FALSE);write.csv(occ,file.path(figdir,"Figure4_source_panel_b.csv"),row.names=FALSE);write.csv(pc_long,file.path(figdir,"Figure4_source_panel_c.csv"),row.names=FALSE);write.csv(gcd,file.path(figdir,"Figure4_source_panel_d.csv"),row.names=FALSE)
cat("NATURE_FIGURE_R_ONLY=PASS\n")
