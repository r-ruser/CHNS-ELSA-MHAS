options(stringsAsFactors=FALSE, cli.unicode=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork);library(svglite);library(ragg)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="C2")stop("Run from C2")
indir<-file.path(project,"results","formal");figdir<-file.path(project,"figures");dir.create(figdir,showWarnings=FALSE)
x<-read.csv(file.path(indir,"main_secondary_results.csv"))
keep_methods<-c("Unadjusted","Minimally adjusted","Fully adjusted","Fully adjusted; community cluster bootstrap","Overlap weighted")
method_labels<-c(
 "Unadjusted"="Unadjusted",
 "Minimally adjusted"="Minimally adjusted",
 "Fully adjusted"="Fully adjusted (sandwich)",
 "Fully adjusted; community cluster bootstrap"="Fully adjusted (bootstrap)",
 "Overlap weighted"="Overlap weighted"
)
src<-bind_rows(
 x%>%filter(outcome=="HD_z",analysis%in%keep_methods)%>%mutate(panel="a",outcome_label="Homeostatic dysregulation",unit="Difference (SD)"),
 x%>%filter(outcome=="KDM_BAA",analysis%in%keep_methods)%>%mutate(panel="b",outcome_label="KDM biological age deviation",unit="Cross-fitted BAA residual difference (years)")
) %>% mutate(
 method=factor(unname(method_labels[analysis]),levels=rev(unname(method_labels[keep_methods]))),
 inference=case_when(
   analysis=="Fully adjusted"~"Primary sandwich",
   analysis=="Fully adjusted; community cluster bootstrap"~"Bootstrap robustness",
   TRUE~"Supporting"
 ),
 estimate_label=sprintf("%+.2f [%+.2f, %+.2f]",estimate,lower,upper))
write.csv(src,file.path(indir,"figure2_source_data.csv"),row.names=FALSE)

pal<-c("Primary sandwich"="#2B6F9C","Bootstrap robustness"="#C77C2B","Supporting"="#8A949E")
theme_nature<-theme_classic(base_size=7,base_family="Arial")+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
 axis.line.x=element_line(linewidth=.35),axis.ticks.x=element_line(linewidth=.35),axis.text=element_text(colour="#222222"),
 axis.title.y=element_blank(),plot.title=element_text(size=8,face="bold",margin=margin(b=5)),plot.subtitle=element_text(size=6.5,colour="#555555"),
 legend.position="none",plot.margin=margin(6,6,5,5))
make_panel<-function(d,title,xlab){
 xr<-range(c(d$lower,d$upper,0));pad<-diff(xr)*.42
 ggplot(d,aes(y=method,x=estimate,colour=inference))+geom_vline(xintercept=0,linetype="22",linewidth=.35,colour="#777777")+
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.12,linewidth=.55)+geom_point(size=2.1)+
  geom_text(aes(x=upper+pad*.08,label=estimate_label),hjust=0,size=2.15,colour="#303030")+
  scale_colour_manual(values=pal)+scale_x_continuous(expand=expansion(mult=c(.07,.47)))+
  labs(title=title,subtitle="Solid-to-clean-only minus persistent any-solid",x=xlab)+theme_nature
}
p1<-make_panel(src%>%filter(panel=="a"),"Homeostatic dysregulation","Adjusted mean difference (SD)")
p2<-make_panel(src%>%filter(panel=="b"),"KDM biological age deviation","Cross-fitted BAA residual difference (years)")
fig<-(p1|p2)+plot_annotation(tag_levels="a",caption=paste(
 "Points are mean differences; bars are 95% CIs. Blue: prespecified fully adjusted community-sandwich CI.",
 "Orange: percentile CI from 1,000 community-cluster bootstrap replicates; constructed outcomes were fixed.",
 sprintf("Analysis denominator: n = %s.",format(unique(src$n)[1],big.mark=",")),sep="\n")) &
 theme(plot.tag=element_text(size=8,face="bold"),plot.caption=element_text(size=6,hjust=0,colour="#444444",margin=margin(t=7),lineheight=.95))
w<-183/25.4;h<-100/25.4;base<-file.path(figdir,"Figure2_fuel_trajectory_biological_age")
svglite::svglite(paste0(base,".svg"),width=w,height=h);print(fig);dev.off()
grDevices::cairo_pdf(paste0(base,".pdf"),width=w,height=h,family="Arial");print(fig);dev.off()
ragg::agg_tiff(paste0(base,".tiff"),width=w,height=h,units="in",res=600,compression="lzw");print(fig);dev.off()
ragg::agg_png(paste0(base,".png"),width=w,height=h,units="in",res=600);print(fig);dev.off()
info<-file.info(paste0(base,c(".svg",".pdf",".tiff",".png")))
stopifnot(all(info$size>1000),nrow(src)>0,all(src$n==unique(src$n)[1]),all(src$lower<=src$estimate&src$estimate<=src$upper))
stopifnot(sum(src$analysis=="Fully adjusted; community cluster bootstrap")==2)
cat("PASS: Figure 2 exported in SVG/PDF/TIFF600/PNG; source rows=",nrow(src),"; bootstrap rows=2\n",sep="")
