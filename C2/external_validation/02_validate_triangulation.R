options(stringsAsFactors=FALSE)
stopifnot(basename(normalizePath(getwd(),winslash="/"))=="external_validation")
fail <- function(x) stop(x,call.=FALSE)
near <- function(a,b,tol=1e-10) all(abs(a-b)<tol,na.rm=TRUE)

required <- c(
  "protocol/triangulation_protocol_v1.0_frozen.md","protocol/triangulation_SAP_v1.0_frozen.md",
  "study_report.md","results/chns_frozen_shared7_parameters.rds",
  "results/elsa_transport_performance.csv","results/elsa_bootstrap_replicates.csv",
  "results/elsa_bootstrap_summary.csv","results/mhas_trajectory_counts.csv",
  "results/triangulation_summary.csv","figures/Figure3_C2_triangular_validation.svg",
  "figures/Figure3_C2_triangular_validation.pdf","figures/Figure3_C2_triangular_validation.tiff",
  "figures/Figure3_C2_triangular_validation.png","figures/figure_legend.md")
if(any(!file.exists(required))) fail(paste("Missing:",paste(required[!file.exists(required)],collapse=", ")))

logtxt <- paste(readLines("logs/01_run_triangulation.log",warn=FALSE,encoding="UTF-8"),collapse="\n")
if(!grepl("TRIANGULATION_RUN=PASS",logtxt,fixed=TRUE)) fail("Run log lacks PASS")
if(grepl("Warning|Error|Execution halted",logtxt,ignore.case=TRUE)) fail("Run log contains warning/error")

p <- read.csv("results/elsa_transport_performance.csv",check.names=FALSE)
if(!identical(p$n,c(1468L,2597L))) fail("Unexpected ELSA denominators")
if(p$pure_transport_pass[1] || p$R2[1]>=.40) fail("ELSA pass status/R2 inconsistent")
b <- read.csv("results/elsa_bootstrap_replicates.csv",check.names=FALSE)
if(nrow(b)!=2000 || any(table(b$version)!=1000)) fail("Bootstrap is not 1000 per version")
if(any(!is.finite(as.matrix(b[,c("MAE","RMSE","R2","calibration_intercept","calibration_slope","BAA_age_correlation","HD7_mean")])))) fail("Non-finite bootstrap metrics")
s <- read.csv("results/elsa_bootstrap_summary.csv",check.names=FALSE)
for(i in seq_len(nrow(s))){
  z <- b[b$version==s$version[i],s$metric[i]]
  q <- as.numeric(quantile(z,c(.025,.975),type=8,names=FALSE))
  if(!near(c(s$lower[i],s$upper[i]),q,1e-9) || !near(s$estimate[i],median(z),1e-9)) fail("Bootstrap summary mismatch")
}

ed <- read.csv("results/elsa_primary_deidentified.csv",check.names=FALSE)
if(nrow(ed)!=1468 || any(grepl("idauniq|hhid|commid|stratum",names(ed),ignore.case=TRUE))) fail("ELSA deidentification check failed")
tc <- read.csv("results/mhas_trajectory_counts.csv",check.names=FALSE)
elig <- unique(tc$eligible_denominator)
if(length(elig)!=1 || elig!=10926 || sum(tc$n[tc$trajectory!="ineligible"])!=elig) fail("MHAS denominator mismatch")
if(any(tc$n[match(c("persistent solid","solid-to-gas"),tc$trajectory)]<200)) fail("MHAS structural rule failed")

srcb <- read.csv("figures/Figure3_source_panel_b.csv",check.names=FALSE)
srcc <- read.csv("figures/Figure3_source_panel_c.csv",check.names=FALSE)
srcd <- read.csv("figures/Figure3_source_panel_d.csv",check.names=FALSE)
if(nrow(srcb)!=1468 || nrow(srcc)!=4 || sum(srcd$n)!=elig) fail("Figure source-data mismatch")

pngdim <- dim(png::readPNG("figures/Figure3_C2_triangular_validation.png"))[1:2]
tifdim <- dim(tiff::readTIFF("figures/Figure3_C2_triangular_validation.tiff"))[1:2]
if(!identical(pngdim,c(1535L,2161L))) fail("PNG dimensions wrong")
if(!identical(tifdim,c(3071L,4323L))) fail("TIFF dimensions wrong")
svg <- paste(readLines("figures/Figure3_C2_triangular_validation.svg",warn=FALSE,encoding="UTF-8"),collapse="\n")
if(length(gregexpr("<text",svg,fixed=TRUE)[[1]])<25) fail("SVG text is not editable")

textfiles <- c(list.files(".",pattern="\\.(md|csv|txt|log)$",recursive=TRUE,full.names=TRUE),
               "figures/Figure3_C2_triangular_validation.svg")
bad <- c("\uFFFD","????","鍏","鏁","鈻","鉁","脳")
for(f in textfiles){
  x <- paste(readLines(f,warn=FALSE,encoding="UTF-8"),collapse="\n")
  if(any(vapply(bad,function(z) grepl(z,x,fixed=TRUE),logical(1)))) fail(paste("Encoding/bad-pattern:",f))
}

manifest_files <- sort(setdiff(list.files(".",recursive=TRUE,full.names=TRUE),"results/final_manifest_md5.csv"))
manifest_files <- manifest_files[file.info(manifest_files)$isdir==FALSE]
manifest <- data.frame(file=sub("^\\./","",manifest_files),md5=unname(tools::md5sum(manifest_files)))
write.csv(manifest,"results/final_manifest_md5.csv",row.names=FALSE,fileEncoding="UTF-8")
cat("C2_TRIANGULATION_VALIDATION=PASS\n")
cat("BOOTSTRAP=1000_PER_ELSA_VERSION\n")
cat("FIGURE_DIMENSIONS=PASS\nENCODING=PASS\n")
