options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dagitty);library(ggdag);library(ggplot2);library(ggraph);library(dplyr);library(svglite)})
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
fd<-file.path(project,"figures","causal");od<-file.path(project,"results","causal");ld<-file.path(project,"logs","causal")
for(z in c(fd,od,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"01_longitudinal_DAG.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);on.exit({sink();close(logcon)},add=TRUE)

dag_main<-dagitty("dag {Age0 -> Fuel_t; Age0 -> HD2009; Sex -> Fuel_t; Sex -> HD2009; Province -> Fuel_t; Province -> HD2009; Baseline_SES -> Fuel_t; Baseline_SES -> SES_t; Baseline_SES -> HD2009; SES_t -> Fuel_t; SES_t -> HD2009; Fuel_t -> Fuel_t1; Fuel_t -> SES_t1; SES_t -> SES_t1; SES_t1 -> Fuel_t1; SES_t1 -> HD2009; Fuel_t1 -> HD2009}")
exposures(dag_main)<-"Fuel_t1";outcomes(dag_main)<-"HD2009"
coordinates(dag_main)<-list(x=c(Age0=0,Sex=0,Province=0,Baseline_SES=0,SES_t=1,Fuel_t=1.3,SES_t1=2.2,Fuel_t1=2.6,HD2009=4),y=c(Age0=3,Sex=2.2,Province=1.4,Baseline_SES=.4,SES_t=.6,Fuel_t=2.3,SES_t1=.8,Fuel_t1=2.2,HD2009=1.8))
dag_supp<-dagitty("dag {Age0 -> Fuel_t; Age0 -> HD2009; Sex -> Fuel_t; Sex -> HD2009; Province -> Fuel_t; Province -> HD2009; Baseline_SES -> Fuel_t; Baseline_SES -> SES_t; Baseline_SES -> HD2009; SES_t -> Fuel_t; SES_t -> HD2009; Fuel_t -> Fuel_t1; Fuel_t -> SES_t1; SES_t -> SES_t1; SES_t1 -> Fuel_t1; SES_t1 -> HD2009; Fuel_t -> Retention_t1; SES_t -> Retention_t1; Retention_t1 -> Biomarker2009; Fuel_t1 -> Biomarker2009; SES_t1 -> Biomarker2009; Biomarker2009 -> HD_observed; Fuel_t1 -> HD2009; HD2009 -> HD_observed}")
exposures(dag_supp)<-"Fuel_t1";outcomes(dag_supp)<-"HD2009"
coordinates(dag_supp)<-list(x=c(Age0=0,Sex=0,Province=0,Baseline_SES=0,SES_t=1,Fuel_t=1.3,SES_t1=2.2,Fuel_t1=2.6,Retention_t1=2.1,Biomarker2009=3.2,HD2009=4,HD_observed=4.7),y=c(Age0=3.2,Sex=2.5,Province=1.8,Baseline_SES=.6,SES_t=.8,Fuel_t=2.5,SES_t1=1,Fuel_t1=2.4,Retention_t1=3.4,Biomarker2009=3.3,HD2009=1.8,HD_observed=2.4))

node_class<-function(n)case_when(grepl("Fuel",n)~"Exposure",grepl("SES",n)~"Time-varying SES",grepl("Retention|Biomarker",n)~"Selection",grepl("HD",n)~"Outcome",TRUE~"Baseline")
plot_dag<-function(g,title){z<-ggdag::tidy_dagitty(g);z$data$class<-node_class(z$data$name);ggdag(z,text=FALSE,use_labels=NULL)+
 geom_dag_node(aes(fill=class),shape=21,size=12,color="white",stroke=.5)+geom_dag_text(aes(label=name),size=2.1,color="#202020")+
 scale_fill_manual(values=c(Baseline="#D9D9D9",Exposure="#5DA5DA","Time-varying SES"="#9FD0C7",Selection="#F2C572",Outcome="#D98585"))+
 theme_dag_blank(base_family="Arial")+theme(legend.position="bottom",legend.title=element_blank(),plot.title=element_text(size=8,face="bold"),plot.margin=margin(5,8,5,8))+labs(title=title)}
save4<-function(p,base,w=183,h=105){svglite(paste0(base,".svg"),w/25.4,h/25.4);print(p);dev.off();cairo_pdf(paste0(base,".pdf"),w/25.4,h/25.4,family="Arial");print(p);dev.off()}
p1<-plot_dag(dag_main,"Longitudinal treatment-confounder feedback");p2<-plot_dag(dag_supp,"Expanded longitudinal DAG with selection")
save4(p1,file.path(fd,"DAG_main"));save4(p2,file.path(fd,"DAG_supplement"),183,120)
writeLines(c("MAIN DAG",as.character(dag_main),"","SUPPLEMENT DAG",as.character(dag_supp)),file.path(od,"dagitty_model.txt"))
adj<-capture.output(print(adjustmentSets(dag_main,type="minimal")));pth<-capture.output(print(paths(dag_main,from="Fuel_t1",to="HD2009")));ici<-capture.output(print(impliedConditionalIndependencies(dag_main)))
writeLines(c("MINIMAL ADJUSTMENT SETS",adj,"","PATHS",pth,"","IMPLIED CONDITIONAL INDEPENDENCIES",ici),file.path(od,"dagitty_queries.txt"))
cat(paste(adj,collapse="\n"),"\nDAG_EXPORT=PASS\n")
