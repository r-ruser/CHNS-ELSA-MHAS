options(stringsAsFactors=FALSE)
stopifnot(basename(normalizePath(getwd(),winslash="/"))=="longitudinal_upgrade")
fail<-function(x)stop(x,call.=FALSE);near<-function(a,b,tol=1e-9)all(abs(a-b)<tol,na.rm=TRUE)
req<-c("protocol/protocol_v1.0_frozen.md","protocol/SAP_v1.0_frozen.md","study_report.md",
"results/transition_intensities.csv","results/occupancy_from_solid.csv","results/dynamic_association_results.csv",
"results/community_bootstrap_replicates.csv","results/community_bootstrap_summary.csv","results/gam_diagnostics.csv",
"results/dynamic_sensitivity_results.csv","data/dynamic_analysis_deidentified.rds",
"figures/Figure4_C2_dynamic_fuel_biological_age.svg","figures/Figure4_C2_dynamic_fuel_biological_age.pdf",
"figures/Figure4_C2_dynamic_fuel_biological_age.tiff","figures/Figure4_C2_dynamic_fuel_biological_age.png",
"figures/figure_legend.md","figures/qa_notes.md")
if(any(!file.exists(req)))fail(paste("Missing",paste(req[!file.exists(req)],collapse=", ")))
lg<-paste(readLines("logs/01_multistate_dynamic_analysis.log",warn=FALSE,encoding="UTF-8"),collapse="\n")
fg<-paste(readLines("logs/03_nature_figure.log",warn=FALSE,encoding="UTF-8"),collapse="\n")
if(!grepl("MULTISTATE_DYNAMIC_ANALYSIS=PASS",lg,fixed=TRUE)||!grepl("NATURE_FIGURE_R_ONLY=PASS",fg,fixed=TRUE))fail("PASS marker missing")
if(grepl("Error|Execution halted",lg,ignore.case=TRUE)||grepl("Warning|Error",fg,ignore.case=TRUE))fail("Final logs contain warning/error")
q<-read.csv("results/transition_intensities.csv");if(nrow(q)!=6||any(q$estimate<0)||any(q$lower>q$estimate|q$upper<q$estimate))fail("Q output invalid")
occ<-read.csv("results/occupancy_from_solid.csv");chk<-aggregate(probability~year,occ,sum);if(max(abs(chk$probability-1))>1e-8)fail("Occupancy probabilities do not sum to one")
ar<-read.csv("results/dynamic_association_results.csv");if(nrow(ar)!=6||any(ar$n!=5149)||any(ar$lower>ar$estimate|ar$upper<ar$estimate))fail("Association output invalid")
b<-read.csv("results/community_bootstrap_replicates.csv");if(nrow(b)!=6000||any(table(b$replicate)!=6)||length(unique(b$replicate))!=1000||any(!is.finite(b$estimate)))fail("Bootstrap invalid")
bs<-read.csv("results/community_bootstrap_summary.csv");for(i in seq_len(nrow(bs))){z<-b$estimate[b$outcome==bs$outcome[i]&b$term==bs$term[i]];qq<-quantile(z,c(.025,.975),type=8,names=FALSE);if(!near(qq,c(bs$lower[i],bs$upper[i])))fail("Bootstrap summary mismatch")}
d<-readRDS("data/dynamic_analysis_deidentified.rds");if(nrow(d)!=5149||any(c("Idind","hhid","commid","stratum","Idind_key")%in%names(d)))fail("Deidentification invalid")
if(nrow(read.csv("figures/Figure4_source_panel_a.csv"))!=6||nrow(read.csv("figures/Figure4_source_panel_b.csv"))!=nrow(occ)||nrow(read.csv("figures/Figure4_source_panel_c.csv"))!=12||nrow(read.csv("figures/Figure4_source_panel_d.csv"))!=202)fail("Figure source data invalid")
pd<-dim(png::readPNG("figures/Figure4_C2_dynamic_fuel_biological_age.png"))[1:2];td<-dim(tiff::readTIFF("figures/Figure4_C2_dynamic_fuel_biological_age.tiff"))[1:2]
if(!identical(pd,c(1772L,2161L))||!identical(td,c(3543L,4323L)))fail("Figure dimensions invalid")
# require nonwhite raster content
im<-png::readPNG("figures/Figure4_C2_dynamic_fuel_biological_age.png");if(mean(im[,,1:3]<.98)<.01)fail("Figure appears blank")
files<-list.files(".",pattern="[.](md|csv|txt|log)$",recursive=TRUE,full.names=TRUE);bad<-c("\uFFFD","????","鍏","鏁","鈻","鉁","脳")
for(f in files){x<-paste(readLines(f,warn=FALSE,encoding="UTF-8"),collapse="\n");if(any(vapply(bad,function(z)grepl(z,x,fixed=TRUE),logical(1))))fail(paste("Encoding",f))}
mf<-sort(setdiff(list.files(".",recursive=TRUE,full.names=TRUE),"results/final_manifest_md5.csv"));mf<-mf[!file.info(mf)$isdir]
write.csv(data.frame(file=mf,md5=unname(tools::md5sum(mf))),"results/final_manifest_md5.csv",row.names=FALSE,fileEncoding="UTF-8")
cat("C2_LONGITUDINAL_UPGRADE_VALIDATION=PASS\nMSM_N=5149\nBOOTSTRAP=1000\nFIGURE_R_ONLY=PASS\nENCODING=PASS\n")
