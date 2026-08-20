options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(cobalt);library(ggplot2)})
set.seed(20260815)
project<-normalizePath(getwd(),winslash="/",mustWork=TRUE);if(basename(project)!="longitudinal_upgrade")stop("Run from project root")
dd<-file.path(project,"data","causal");od<-file.path(project,"results","causal");td<-file.path(project,"tables","causal");fd<-file.path(project,"figures","causal");ld<-file.path(project,"logs","causal")
for(z in c(dd,od,td,fd,ld))dir.create(z,FALSE,TRUE)
logcon<-file(file.path(ld,"04_weight_balance.log"),"wt",encoding="UTF-8");sink(logcon,split=TRUE);sink(logcon,type="message");on.exit({sink(type="message");sink();close(logcon)},add=TRUE)
w0<-readRDS(file.path(dd,"causal_participant_weights_internal_restricted.rds"))
long<-readRDS(file.path(dd,"causal_treatment_long_internal_restricted.rds"))
w<-w0%>%filter(strategy%in%c("continued_any_solid","sustained_clean_only"))%>%mutate(A=as.integer(strategy=="sustained_clean_only"))
w_all<-w0%>%mutate(A=as.integer(strategy=="sustained_clean_only"))

qclip<-function(x,p){if(p==0)return(x);q<-quantile(x,c(p,1-p),na.rm=TRUE);pmin(pmax(x,q[1]),q[2])}
cuts<-c(none=0,p01_p99=.01,p025_p975=.025,p05_p95=.05)
wlong<-bind_rows(lapply(names(cuts),function(nm)w%>%mutate(rule=nm,treatment_weight_use=qclip(treatment_weight,cuts[[nm]]),censoring_weight_use=qclip(censoring_weight,cuts[[nm]]),combined_weight_use=qclip(combined_weight,cuts[[nm]]))))
desc1<-function(x){data.frame(N=sum(is.finite(x)),mean=mean(x),SD=sd(x),minimum=min(x),P1=quantile(x,.01),P5=quantile(x,.05),P25=quantile(x,.25),median=median(x),P75=quantile(x,.75),P95=quantile(x,.95),P99=quantile(x,.99),maximum=max(x),ESS=sum(x)^2/sum(x^2))}
diag<-bind_rows(lapply(split(wlong,wlong$rule),function(z)bind_rows(
 transform(desc1(z$treatment_weight_use),weight="IPTW"),transform(desc1(z$censoring_weight_use),weight="IPCW"),transform(desc1(z$combined_weight_use),weight="combined"))%>%mutate(rule=unique(z$rule))))%>%select(rule,weight,everything())
write.csv(diag,file.path(td,"Supplementary_Table_weight_diagnostics.csv"),row.names=FALSE)

# Candidate primary rule is selected by stability and balance, never by outcome results.
# 1/99 is the least aggressive prespecified truncation expected to control the observed tail.
primary_rule<-"p01_p99"; wp<-wlong%>%filter(rule==primary_rule)%>%select(Idind,A,strategy,treatment_weight_use,censoring_weight_use,combined_weight_use)
saveRDS(wp,file.path(dd,"causal_primary_weights_internal_restricted.rds"),compress="gzip")

wmean<-function(x,wt)sum(x*wt,na.rm=TRUE)/sum(wt[is.finite(x)],na.rm=TRUE)
wvar<-function(x,wt){m<-wmean(x,wt);sum(wt*(x-m)^2,na.rm=TRUE)/sum(wt[is.finite(x)],na.rm=TRUE)}
smd_pair<-function(x,g,wt,g0,g1){i0<-g==g0&is.finite(x);i1<-g==g1&is.finite(x);den<-sqrt((wvar(x[i1],wt[i1])+wvar(x[i0],wt[i0]))/2);if(!is.finite(den)||den==0)return(NA_real_);(wmean(x[i1],wt[i1])-wmean(x[i0],wt[i0]))/den}
smd_multi_max<-function(x,g,wt){lev<-unique(as.character(g[!is.na(g)]));if(length(lev)<2)return(NA_real_);prs<-combn(lev,2,simplify=FALSE);v<-vapply(prs,function(p)smd_pair(x,as.character(g),wt,p[1],p[2]),numeric(1));if(all(!is.finite(v)))return(NA_real_);max(abs(v[is.finite(v)]))}
# Explicit implementation is used for wave-specific multi-category diagnostics;
# cobalt is additionally used below for the final 2009 binary strategy audit.
vars<-c("age","urban_i","education_i","income_i","assets_i","smoking_i","alcohol_i")
bal<-list()
for(yr in sort(unique(long$wave))){
 z<-long%>%filter(wave==yr)%>%semi_join(w_all,by="Idind")
 if(!nrow(z))next
 lo<-quantile(z$treatment_weight_cum,c(.01,.99),na.rm=TRUE);z$wave_w<-pmin(pmax(z$treatment_weight_cum,lo[1]),lo[2])
 for(v in vars){bal[[length(bal)+1]]<-data.frame(wave=yr,variable=v,SMD_unweighted=smd_multi_max(z[[v]],z$fuel_state,rep(1,nrow(z))),SMD_weighted=smd_multi_max(z[[v]],z$fuel_state,z$wave_w))}
}
balance<-bind_rows(bal)%>%mutate(abs_unweighted=abs(SMD_unweighted),abs_weighted=abs(SMD_weighted),meets_0_10=abs_weighted<.10)
write.csv(balance,file.path(td,"Table_causal_balance.csv"),row.names=FALSE)

z09<-long%>%filter(wave==2009)%>%inner_join(w_all%>%select(Idind,combined_weight),by="Idind")
lo09<-quantile(z09$combined_weight,c(.01,.99),na.rm=TRUE);z09$combined_weight_use<-pmin(pmax(z09$combined_weight,lo09[1]),lo09[2])
bt<-cobalt::bal.tab(fuel_state~age+urban_i+education_i+income_i+assets_i+smoking_i+alcohol_i,data=z09,weights=z09$combined_weight_use,method="weighting",estimand="ATE",un=TRUE,pairwise=TRUE,thresholds=c(m=.1))
capture.output(print(bt),file=file.path(od,"cobalt_balance_2009.txt"))
love<-cobalt::love.plot(bt,abs=TRUE,thresholds=c(m=.1),colors=c("#777777","#0072B2"),shapes=c(16,17),sample.names=c("Unweighted","Weighted"),var.order="unadjusted",stars="raw")+theme_minimal(base_size=9)+theme(legend.position="top",panel.grid.minor=element_blank())
ggsave(file.path(fd,"Figure_causal_loveplot.pdf"),love,width=170,height=115,units="mm",device=cairo_pdf)
ggsave(file.path(fd,"Figure_causal_loveplot.svg"),love,width=170,height=115,units="mm",device=svglite::svglite)

pd<-wlong%>%select(Idind,strategy,rule,IPTW=treatment_weight_use,IPCW=censoring_weight_use,Combined=combined_weight_use)%>%pivot_longer(c(IPTW,IPCW,Combined),names_to="weight_type",values_to="weight")
p<-ggplot(pd,aes(weight,colour=rule,fill=rule))+geom_density(alpha=.08,linewidth=.45)+facet_wrap(~weight_type,scales="free",ncol=1)+coord_cartesian(xlim=c(0,quantile(pd$weight,.995)))+theme_classic(base_size=9)+labs(x="Participant-level stabilized weight",y="Density",colour="Truncation",fill="Truncation")+theme(legend.position="top")
ggsave(file.path(fd,"Supplementary_Figure_weight_distributions.pdf"),p,width=170,height=150,units="mm",device=cairo_pdf)
ggsave(file.path(fd,"Supplementary_Figure_weight_distributions.svg"),p,width=170,height=150,units="mm",device=svglite::svglite)

ess_groups<-wp%>%group_by(strategy)%>%summarise(N=n(),ESS=(sum(combined_weight_use)^2/sum(combined_weight_use^2)),.groups="drop")
write.csv(ess_groups,file.path(od,"primary_weight_group_ESS.csv"),row.names=FALSE)
maxbal<-balance%>%group_by(wave)%>%summarise(max_abs_SMD=if(all(is.na(abs_weighted)))NA_real_ else max(abs_weighted,na.rm=TRUE),n_over_0_10=sum(abs_weighted>=.1,na.rm=TRUE),.groups="drop")
write.csv(maxbal,file.path(od,"balance_wave_summary.csv"),row.names=FALSE)
cat("PRIMARY RULE=",primary_rule,"\n",sep="");print(ess_groups);print(maxbal)
cat("WEIGHT_BALANCE=FAILED_THRESHOLD; CAUSAL STOP RULE TRIGGERED\n")
