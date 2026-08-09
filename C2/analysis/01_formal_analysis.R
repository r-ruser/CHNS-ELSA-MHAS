options(stringsAsFactors = FALSE, warn = 1, cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(splines)
  library(sandwich); library(lmtest); library(broom)
})

project <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (basename(project) != "C2") stop("Run with working directory set to C2")
workspace <- dirname(project)
dir.create(file.path(project, "data"), showWarnings = FALSE)
dir.create(file.path(project, "results", "formal"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project, "logs"), showWarnings = FALSE)
outdir <- file.path(project, "results", "formal")
logfile <- file.path(project, "logs", "03_formal_analysis.log")
logcon <- file(logfile, open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE); sink(logcon, type = "message")
on.exit({sink(type = "message"); sink(); close(logcon)}, add = TRUE)
cat("C2 formal analysis started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

find_one <- function(pattern) {
  z <- list.files(workspace, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  z <- z[grepl("CHNS-STATA", z, fixed = TRUE)]
  if (length(z) != 1L) stop("Expected one CHNS-STATA match for ", pattern, "; found ", length(z))
  normalizePath(z, winslash = "/", mustWork = TRUE)
}
files <- c(
  asset = find_one("^asset_12\\.dta$"), biomarker = find_one("^biomarker_09\\.dta$"),
  surveys = find_one("^surveys_pub_12\\.dta$"), master = find_one("^mast_pub_12\\.dta$"),
  pexam = find_one("^pexam_pub_12\\.dta$"), education = find_one("^educ_12\\.dta$"),
  hh_income = find_one("^hhinc_10\\.dta$")
)
write.csv(data.frame(dataset = names(files), path = unname(files)), file.path(outdir, "formal_input_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")

asset <- read_dta(files["asset"]) %>% select(hhid, wave, L8_1, L8_2, commid, T1,
  any_of(c("L23", "L27", "L31", "L105", "L110", "L115", "L120", "L140E")))
surveys <- read_dta(files["surveys"]) %>% select(Idind, hhid, wave, age, commid, t1, stratum, urban)
bio <- read_dta(files["biomarker"])
master <- read_dta(files["master"]) %>% select(Idind, gender)
pexam <- read_dta(files["pexam"]) %>% filter(wave == 2009) %>%
  select(IDind, U25, U40, U22, U24A, SYSTOL1, SYSTOL2, SYSTOL3, DIASTOL1, DIASTOL2, DIASTOL3, height, weight, U10)
educ <- read_dta(files["education"]) %>% filter(wave == 2009) %>% select(IDind, A11, A12)
income <- read_dta(files["hh_income"]) %>% filter(wave == 2009) %>% select(hhid, hhincpc_cpi)

stopifnot(nrow(asset %>% filter(wave <= 2009) %>% count(hhid, wave) %>% filter(n > 1)) == 0)
stopifnot(nrow(surveys %>% count(Idind, wave) %>% filter(n > 1)) == 0)
stopifnot(nrow(bio %>% count(IDind) %>% filter(n > 1)) == 0)

class_primary <- function(x) case_when(x %in% c(1, 6, 7) ~ "solid", x %in% c(2, 4, 5) ~ "clean", TRUE ~ "other")
fuel_long <- surveys %>% filter(wave <= 2009) %>%
  left_join(asset %>% filter(wave <= 2009) %>% select(hhid, wave, L8_1, L8_2), by = c("hhid", "wave")) %>%
  mutate(primary_state = class_primary(L8_1), secondary_state = case_when(
    is.na(L8_2) | L8_2 == 0 ~ "none", L8_2 %in% c(1, 6, 7) ~ "solid",
    L8_2 %in% c(2, 4, 5) ~ "clean", TRUE ~ "other"),
    stacking = primary_state %in% c("solid", "clean") & secondary_state %in% c("solid", "clean") & primary_state != secondary_state,
    effective_state = case_when(primary_state == "solid" | secondary_state == "solid" ~ "any_solid",
      primary_state == "clean" & secondary_state %in% c("none", "clean") ~ "clean_only", TRUE ~ NA_character_))

classify_sequence <- function(state, wave, solid_label, clean_label) {
  o <- order(wave); state <- state[o]; state <- state[!is.na(state)]
  if (!length(state)) return("no_valid_wave")
  if (length(state) < 2) return("fewer_than_2_valid_waves")
  if (all(state == solid_label)) return("persistent_any_solid")
  if (all(state == clean_label)) return("persistent_clean_only")
  tr <- sum(state[-1] != state[-length(state)])
  if (state[1] == solid_label && tail(state, 1) == clean_label && tr == 1) return("solid_to_clean_only")
  if (state[1] == clean_label && tail(state, 1) == solid_label && tr == 1) return("clean_to_solid")
  "complex_switching"
}
first_valid <- function(state,wave){o<-order(wave);z<-state[o];z<-z[!is.na(z)];if(length(z))z[1]else NA_character_}
last_valid <- function(state,wave){o<-order(wave);z<-state[o];z<-z[!is.na(z)];if(length(z))tail(z,1)else NA_character_}

traj_effective <- fuel_long %>% group_by(Idind) %>% summarise(
  n_prior_effective = sum(wave < 2009 & !is.na(effective_state)),
  known_effective_2009 = any(wave == 2009 & !is.na(effective_state)),
  stacking_ever = any(stacking, na.rm = TRUE),
  first_effective_state=first_valid(effective_state,wave),last_effective_state=last_valid(effective_state,wave),
  trajectory = classify_sequence(effective_state, wave, "any_solid", "clean_only"), .groups = "drop")
traj_primary <- fuel_long %>% mutate(pstate = ifelse(primary_state == "solid", "any_solid", ifelse(primary_state == "clean", "clean_only", NA))) %>%
  group_by(Idind) %>% summarise(n_prior_primary = sum(wave < 2009 & !is.na(pstate)), known_primary_2009 = any(wave == 2009 & !is.na(pstate)),
    first_primary_state=first_valid(pstate,wave),last_primary_state=last_valid(pstate,wave),
    trajectory_primary_only = classify_sequence(pstate, wave, "any_solid", "clean_only"), .groups = "drop")

asset_vars <- intersect(c("L23", "L27", "L31", "L105", "L110", "L115", "L120", "L140E"), names(asset))
a09 <- asset %>% filter(wave == 2009) %>% mutate(across(all_of(asset_vars), ~ifelse(.x %in% c(0, 1), .x, NA_real_)))
a09$asset_count <- apply(as.data.frame(a09[, asset_vars, drop = FALSE]), 1, function(z) if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE))
a09 <- a09 %>% select(hhid, asset_count)
s09 <- surveys %>% filter(wave == 2009) %>% select(Idind, hhid, age, commid, province = t1, stratum, urban)
p09 <- pexam %>% transmute(Idind = IDind, smoking = ifelse(U25 %in% c(0,1), U25, NA_real_), alcohol = ifelse(U40 %in% c(0,1), U40, NA_real_),
  hypertension = ifelse(U22 %in% c(0,1), U22, NA_real_), diabetes = ifelse(U24A %in% c(0,1), U24A, NA_real_),
  sbp = rowMeans(cbind(SYSTOL1, SYSTOL2, SYSTOL3), na.rm = TRUE), dbp = rowMeans(cbind(DIASTOL1, DIASTOL2, DIASTOL3), na.rm = TRUE),
  bmi = weight/(height/100)^2, waist = U10)
p09$sbp[!is.finite(p09$sbp)] <- NA; p09$dbp[!is.finite(p09$dbp)] <- NA
e09 <- educ %>% transmute(Idind = IDind, education_year_code = ifelse(A11 >= 0, A11, NA_real_), education_level_code = ifelse(A12 %in% 0:6, A12, NA_real_))

dat <- bio %>% rename(Idind = IDind) %>% left_join(s09, by = "Idind") %>% left_join(master, by = "Idind") %>%
  left_join(p09, by = "Idind") %>% left_join(e09, by = "Idind") %>% left_join(income, by = "hhid") %>%
  left_join(a09, by = "hhid") %>% left_join(traj_effective, by = "Idind") %>% left_join(traj_primary, by = "Idind")
stopifnot(nrow(dat) == nrow(bio))

raw_bm <- c("HS_CRP", "alb", "cre", "glucose", "tg", "HDL_C", "wbc", "hgb", "sbp")
trans_names <- c("ln_hscrp", "albumin", "ln_creatinine", "ln_glucose", "ln_triglyceride", "hdl", "ln_wbc", "hemoglobin", "sbp")
transform_bm <- function(d) data.frame(
  ln_hscrp = log1p(pmax(d$HS_CRP, 0)), albumin = d$alb, ln_creatinine = log(pmax(d$cre, 1e-6)),
  ln_glucose = log(pmax(d$glucose, 1e-6)), ln_triglyceride = log(pmax(d$tg, 1e-6)), hdl = d$HDL_C,
  ln_wbc = log(pmax(d$wbc, 1e-6)), hemoglobin = d$hgb, sbp = d$sbp)
Xraw <- transform_bm(dat)
adult_complete <- !is.na(dat$age) & dat$age >= 18 & !is.na(dat$gender) & complete.cases(Xraw)

winsor_limits <- function(x, p = c(.005, .995)) quantile(x, p, na.rm = TRUE, names = FALSE, type = 8)
lims <- lapply(Xraw[adult_complete, , drop = FALSE], winsor_limits)
X <- Xraw
for (j in names(X)) X[[j]] <- pmin(pmax(X[[j]], lims[[j]][1]), lims[[j]][2])

# Primary HD: pre-specified transformations; age 20-39 reference; median center;
# covariance regularized by a very small diagonal ridge; standardized log distance.
ref <- adult_complete & dat$age >= 20 & dat$age <= 39
if (sum(ref) < 500) stop("HD reference sample unexpectedly small")
center <- sapply(X[ref, , drop = FALSE], median, na.rm = TRUE)
S <- cov(X[ref, , drop = FALSE]); ridge <- median(diag(S)) * 1e-6; Sinv <- solve(S + diag(ridge, ncol(S)))
md2 <- rep(NA_real_, nrow(dat)); dx <- sweep(as.matrix(X[adult_complete, , drop = FALSE]), 2, center)
md2[adult_complete] <- rowSums((dx %*% Sinv) * dx)
hd_raw <- log1p(sqrt(pmax(md2, 0)))
dat$HD_z <- as.numeric((hd_raw - mean(hd_raw[adult_complete], na.rm = TRUE))/sd(hd_raw[adult_complete], na.rm = TRUE))

# Alternative HD uses the full adult reference and 1/99% winsorization.
Xalt <- Xraw; lims_alt <- lapply(Xraw[adult_complete, , drop = FALSE], winsor_limits, p = c(.01,.99))
for (j in names(Xalt)) Xalt[[j]] <- pmin(pmax(Xalt[[j]], lims_alt[[j]][1]), lims_alt[[j]][2])
c2 <- sapply(Xalt[adult_complete, , drop = FALSE], median); S2 <- cov(Xalt[adult_complete, , drop = FALSE]); inv2 <- solve(S2 + diag(median(diag(S2))*1e-6, ncol(S2)))
dx2 <- sweep(as.matrix(Xalt[adult_complete, , drop = FALSE]), 2, c2)
hd2 <- rep(NA_real_, nrow(dat)); hd2[adult_complete] <- log1p(sqrt(pmax(rowSums((dx2 %*% inv2)*dx2), 0)))
dat$HD_alt_z <- as.numeric((hd2-mean(hd2[adult_complete],na.rm=TRUE))/sd(hd2[adult_complete],na.rm=TRUE))

# Five-fold cross-fitted, sex-stratified KDM. Every estimated limit, biomarker-age
# regression, KDM variance and BAA residual model is fitted without the held-out fold.
set.seed(20260809)
pool <- which(adult_complete)
age_band <- cut(dat$age, breaks = c(18,30,40,50,60,70,80,Inf), right = FALSE)
fold <- rep(NA_integer_, nrow(dat))
strata <- interaction(dat$gender, age_band, drop = TRUE)
for (ss in levels(strata)) { ii <- pool[strata[pool] == ss]; fold[ii] <- sample(rep(1:5, length.out = length(ii))) }
kdm_ba <- kdm_baa <- rep(NA_real_, nrow(dat)); kdm_diag <- list()
for (f in 1:5) for (sx in sort(unique(dat$gender[pool]))) {
  tr <- pool[fold[pool] != f & dat$gender[pool] == sx]; te <- pool[fold[pool] == f & dat$gender[pool] == sx]
  if (length(te) == 0 || length(tr) < 200) next
  Xtr <- Xraw[tr,,drop=FALSE]; Xte <- Xraw[te,,drop=FALSE]
  flims <- lapply(Xtr, winsor_limits)
  for (j in names(Xtr)) { Xtr[[j]] <- pmin(pmax(Xtr[[j]], flims[[j]][1]), flims[[j]][2]); Xte[[j]] <- pmin(pmax(Xte[[j]], flims[[j]][1]), flims[[j]][2]) }
  q <- k <- se <- numeric(ncol(Xtr))
  for (j in seq_along(Xtr)) { fitj <- lm(Xtr[[j]] ~ dat$age[tr]); q[j] <- coef(fitj)[1]; k[j] <- coef(fitj)[2]; se[j] <- sigma(fitj) }
  den <- sum(k^2/se^2)
  bae_tr <- as.numeric((as.matrix(sweep(Xtr, 2, q)) %*% (k/se^2))/den)
  bae_te <- as.numeric((as.matrix(sweep(Xte, 2, q)) %*% (k/se^2))/den)
  sba2 <- var(bae_tr - dat$age[tr]); if (!is.finite(sba2) || sba2 < 1e-6) sba2 <- 1e-6
  ba_tr <- (bae_tr*den + dat$age[tr]/sba2)/(den + 1/sba2)
  ba_te <- (bae_te*den + dat$age[te]/sba2)/(den + 1/sba2)
  accfit <- lm(ba_tr ~ dat$age[tr])
  baa_te <- ba_te - (coef(accfit)[1] + coef(accfit)[2]*dat$age[te])
  kdm_ba[te] <- ba_te; kdm_baa[te] <- baa_te
  calfit<-lm(dat$age[te]~ba_te)
  kdm_diag[[length(kdm_diag)+1]] <- data.frame(fold=f, sex=sx, n_train=length(tr), n_test=length(te), sba2=sba2,
    age_rmse=sqrt(mean((ba_te-dat$age[te])^2)),age_mae=mean(abs(ba_te-dat$age[te])),age_correlation=cor(ba_te,dat$age[te]),
    age_r_squared=summary(calfit)$r.squared,calibration_intercept=coef(calfit)[1],calibration_slope=coef(calfit)[2],baa_age_correlation=cor(baa_te,dat$age[te]))
}
dat$KDM_BA <- kdm_ba; dat$KDM_BAA <- kdm_baa
dat$KDM_BAA_z <- as.numeric((kdm_baa-mean(kdm_baa[pool],na.rm=TRUE))/sd(kdm_baa[pool],na.rm=TRUE))
dat$KDM_BAA_direct_z <- as.numeric(((kdm_ba-dat$age)-mean((kdm_ba-dat$age)[pool],na.rm=TRUE))/sd((kdm_ba-dat$age)[pool],na.rm=TRUE))
kdmd<-bind_rows(kdm_diag);write.csv(kdmd, file.path(outdir,"kdm_crossfit_diagnostics.csv"), row.names=FALSE)
kdm_thresholds<-data.frame(
 diagnostic=c("Minimum training n","Maximum MAE (years)","Minimum R-squared","Calibration slope range","Maximum absolute BAA-age correlation"),
 observed=c(min(kdmd$n_train),max(kdmd$age_mae),min(kdmd$age_r_squared),sprintf("%.3f to %.3f",min(kdmd$calibration_slope),max(kdmd$calibration_slope)),max(abs(kdmd$baa_age_correlation))),
 threshold=c(">=200","<=15",">=0.50","0.50 to 2.00","<0.10"),
 pass=c(min(kdmd$n_train)>=200,max(kdmd$age_mae)<=15,min(kdmd$age_r_squared)>=.50,min(kdmd$calibration_slope)>=.5&max(kdmd$calibration_slope)<=2,max(abs(kdmd$baa_age_correlation))<.10))
write.csv(kdm_thresholds,file.path(outdir,"kdm_diagnostic_thresholds.csv"),row.names=FALSE)

dat <- dat %>% mutate(
  sex = factor(gender, levels=c(1,2), labels=c("Male","Female")), province_f=factor(ifelse(is.na(province),"Missing",as.character(province))),
  urban_f=factor(ifelse(is.na(urban),"Missing",ifelse(urban==1,"Urban/town","Rural/suburban"))),
  education = case_when(education_level_code %in% c(0,1) ~ "Primary or less", education_level_code==2 ~ "Lower middle",
    education_level_code %in% c(3,4) ~ "Upper/vocational", education_level_code %in% c(5,6) ~ "College+", TRUE ~ "Missing"),
  smoking_f = factor(case_when(smoking==1 ~ "Yes", smoking==0 ~ "No", TRUE ~ "Missing")),
  alcohol_f = factor(case_when(alcohol==1 ~ "Yes", alcohol==0 ~ "No", TRUE ~ "Missing")),
  trajectory = factor(trajectory, levels=c("persistent_any_solid","solid_to_clean_only","persistent_clean_only","complex_switching","clean_to_solid")),
  income_asinh = asinh(hhincpc_cpi/1000), income_missing=as.integer(is.na(income_asinh)), asset_missing=as.integer(is.na(asset_count)))
dat$income_asinh[is.na(dat$income_asinh)] <- median(dat$income_asinh,na.rm=TRUE)
dat$asset_count[is.na(dat$asset_count)] <- median(dat$asset_count,na.rm=TRUE)

eligible <- with(dat, age>=18 & !is.na(age) & !is.na(sex) & known_effective_2009 & n_prior_effective>=2 & !is.na(HD_z) & !is.na(KDM_BAA_z))
analysis <- dat[eligible,]
analysis$Idind_key <- sprintf("%.0f", analysis$Idind)
saveRDS(analysis, file.path(project,"data","analysis_dataset.rds"), compress="gzip") # legacy internal filename
saveRDS(analysis, file.path(project,"data","analysis_dataset_internal_restricted.rds"), compress="gzip")
analysis_deid<-analysis%>%arrange(Idind)%>%mutate(record_id=sprintf("C2-%05d",row_number()))%>%
  select(record_id,everything(),-any_of(c("Idind","hhid","commid","stratum","Idind_key")))
saveRDS(analysis_deid,file.path(project,"data","analysis_dataset_deidentified.rds"),compress="gzip")
write.csv(analysis_deid,file.path(project,"data","analysis_dataset_deidentified.csv"),row.names=FALSE,fileEncoding="UTF-8")

selection_cohort<-dat[adult_complete,]%>%mutate(selected_formal=eligible[adult_complete])
saveRDS(selection_cohort,file.path(project,"data","selection_cohort_internal_restricted.rds"),compress="gzip")
sel_summary<-bind_rows(lapply(c(FALSE,TRUE),function(z){dd<-selection_cohort[selection_cohort$selected_formal==z,];data.frame(selected_formal=z,n=nrow(dd),age_mean=mean(dd$age),female_pct=100*mean(dd$sex=="Female"),urban_pct=100*mean(dd$urban==1,na.rm=TRUE),income_mean=mean(dd$income_asinh),asset_mean=mean(dd$asset_count),HD_mean=mean(dd$HD_z),KDM_BAA_mean=mean(dd$KDM_BAA))}))
write.csv(sel_summary,file.path(outdir,"selection_characteristics.csv"),row.names=FALSE)

flow <- data.frame(stage=c("2009 biomarker records","Linked age and sex; adult","Complete locked 9-system biomarkers","Known 2009 joint fuel state and >=2 prior effective waves","Formal analysis dataset","Primary two-trajectory contrast"),
 n=c(nrow(dat),sum(!is.na(dat$age)&dat$age>=18&!is.na(dat$gender)),sum(adult_complete),sum(eligible),nrow(analysis),sum(analysis$trajectory %in% c("persistent_any_solid","solid_to_clean_only"))))
write.csv(flow,file.path(outdir,"formal_sample_flow.csv"),row.names=FALSE)

covars <- c("age","sex","province_f","urban_f","education","income_asinh","income_missing","asset_count","asset_missing","smoking_f","alcohol_f")
miss <- data.frame(variable=c(raw_bm,covars,"hypertension","diabetes"), n_missing=sapply(c(raw_bm,covars,"hypertension","diabetes"),function(v)sum(is.na(analysis[[v]]))))
miss$pct_missing <- 100*miss$n_missing/nrow(analysis); write.csv(miss,file.path(outdir,"missingness.csv"),row.names=FALSE)

primary <- analysis %>% filter(trajectory %in% c("persistent_any_solid","solid_to_clean_only")) %>% droplevels()
primary$treated <- as.integer(primary$trajectory=="solid_to_clean_only")
full_rhs <- "trajectory + ns(age,4) + sex + province_f + urban_f + education + income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
min_rhs <- "trajectory + ns(age,4) + sex + province_f + urban_f"

vc_cluster <- function(fit,d,level="community") {
  cl <- if(level=="community") d$commid else d$hhid
  sandwich::vcovCL(fit, cluster=cl, type="HC1", fix=TRUE)
}
std_contrast <- function(outcome,rhs,d,model_name,cluster="community",a="solid_to_clean_only",b="persistent_any_solid") {
  fit <- lm(as.formula(paste(outcome,"~",rhs)),data=d)
  V <- vc_cluster(fit,d,cluster); nd1<-d; nd0<-d; nd1$trajectory<-factor(a,levels=levels(d$trajectory)); nd0$trajectory<-factor(b,levels=levels(d$trajectory))
  tt<-delete.response(terms(fit)); X1<-model.matrix(tt,nd1,contrasts.arg=fit$contrasts,xlev=fit$xlevels); X0<-model.matrix(tt,nd0,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
  cc<-colMeans(X1-X0); bb<-coef(fit); keep<-intersect(intersect(names(cc),names(bb)[!is.na(bb)]),colnames(V)); cc<-cc[keep]; bb<-bb[keep]; V<-V[keep,keep,drop=FALSE]
  est<-sum(cc*bb); se<-sqrt(as.numeric(t(cc)%*%V%*%cc))
  data.frame(outcome=outcome,analysis=model_name,contrast=paste(a,"vs",b),estimate=est,se=se,lower=est-1.96*se,upper=est+1.96*se,p_value=2*pnorm(-abs(est/se)),n=nrow(d),clusters=n_distinct(if(cluster=="community")d$commid else d$hhid),cluster=cluster)
}
models <- bind_rows(
  std_contrast("HD_z","trajectory",primary,"Unadjusted"), std_contrast("HD_z",min_rhs,primary,"Minimally adjusted"), std_contrast("HD_z",full_rhs,primary,"Fully adjusted"),
  std_contrast("KDM_BAA_z","trajectory",primary,"Unadjusted"), std_contrast("KDM_BAA_z",min_rhs,primary,"Minimally adjusted"), std_contrast("KDM_BAA_z",full_rhs,primary,"Fully adjusted"),
  std_contrast("KDM_BAA","trajectory",primary,"Unadjusted"), std_contrast("KDM_BAA",min_rhs,primary,"Minimally adjusted"), std_contrast("KDM_BAA",full_rhs,primary,"Fully adjusted"),
  std_contrast("HD_z",full_rhs,primary,"Fully adjusted; household clustered","household"), std_contrast("KDM_BAA_z",full_rhs,primary,"Fully adjusted; household clustered","household")
)

# Propensity overlap diagnostics and weighted sensitivity analysis.
psfit <- glm(treated ~ ns(age,4)+sex+province_f+urban_f+education+income_asinh+income_missing+asset_count+asset_missing+smoking_f+alcohol_f, family=binomial(), data=primary)
primary$ps <- pmin(pmax(predict(psfit,type="response"),1e-4),1-1e-4)
primary$ow <- ifelse(primary$treated==1,1-primary$ps,primary$ps)
ess <- function(w) sum(w)^2/sum(w^2)
wd <- primary %>% group_by(treated) %>% summarise(n=n(),ps_min=min(ps),ps_q01=quantile(ps,.01),ps_median=median(ps),ps_q99=quantile(ps,.99),ps_max=max(ps),
  weight_mean=mean(ow),weight_sd=sd(ow),weight_p05=quantile(ow,.05),weight_p95=quantile(ow,.95),weight_p99=quantile(ow,.99),weight_max=max(ow),ESS=ess(ow),.groups="drop")
write.csv(wd,file.path(outdir,"propensity_weight_diagnostics.csv"),row.names=FALSE)
write.csv(data.frame(weight_type="overlap",truncation_applied=FALSE,truncation_rule="None; bounded overlap weights used",minimum=min(primary$ow),maximum=max(primary$ow)),file.path(outdir,"weight_truncation_contract.csv"),row.names=FALSE)
overlap <- data.frame(common_support_lower=max(tapply(primary$ps,primary$treated,min)),common_support_upper=min(tapply(primary$ps,primary$treated,max)),
  n_outside_common_support=sum(primary$ps<max(tapply(primary$ps,primary$treated,min))|primary$ps>min(tapply(primary$ps,primary$treated,max))))
write.csv(overlap,file.path(outdir,"positivity_overlap.csv"),row.names=FALSE)

w_smd <- function(x,a,w=NULL) { ok<-!is.na(x)&!is.na(a);x<-x[ok];a<-a[ok];if(is.null(w))w<-rep(1,length(x))else w<-w[ok];
 m1<-weighted.mean(x[a==1],w[a==1]);m0<-weighted.mean(x[a==0],w[a==0]);v1<-weighted.mean((x[a==1]-m1)^2,w[a==1]);v0<-weighted.mean((x[a==0]-m0)^2,w[a==0]);(m1-m0)/sqrt((v1+v0)/2)}
bal_vars <- c("age","urban","income_asinh","asset_count","smoking","alcohol")
balance <- bind_rows(lapply(bal_vars,function(v)data.frame(variable=v,level="continuous/binary",smd_unweighted=w_smd(primary[[v]],primary$treated),smd_overlap=w_smd(primary[[v]],primary$treated,primary$ow))))
for(v in c("sex","education","province_f")) for(lv in levels(factor(primary[[v]]))) balance<-bind_rows(balance,data.frame(variable=v,level=lv,smd_unweighted=w_smd(as.numeric(primary[[v]]==lv),primary$treated),smd_overlap=w_smd(as.numeric(primary[[v]]==lv),primary$treated,primary$ow)))
write.csv(balance,file.path(outdir,"covariate_balance.csv"),row.names=FALSE)

weighted_result <- function(outcome) { fit<-lm(as.formula(paste(outcome,"~treated")),data=primary,weights=ow);V<-vc_cluster(fit,primary,"community");e<-coef(fit)["treated"];s<-sqrt(V["treated","treated"]);data.frame(outcome=outcome,analysis="Overlap weighted",contrast="solid_to_clean_only vs persistent_any_solid",estimate=e,se=s,lower=e-1.96*s,upper=e+1.96*s,p_value=2*pnorm(-abs(e/s)),n=nrow(primary),clusters=n_distinct(primary$commid),cluster="community") }
models<-bind_rows(models,weighted_result("HD_z"),weighted_result("KDM_BAA_z"),weighted_result("KDM_BAA"));write.csv(models,file.path(outdir,"main_secondary_results.csv"),row.names=FALSE)

# Table 1 (formal primary contrast).
qfmt<-function(x)sprintf("%.2f (%.2f)",mean(x,na.rm=TRUE),sd(x,na.rm=TRUE)); cfmt<-function(x)sprintf("%d (%.1f%%)",sum(x,na.rm=TRUE),100*mean(x,na.rm=TRUE))
tab1<-bind_rows(lapply(levels(primary$trajectory),function(g){d<-primary[primary$trajectory==g,];data.frame(group=g,n=nrow(d),age=qfmt(d$age),female=cfmt(d$sex=="Female"),urban=cfmt(d$urban==1),college=cfmt(d$education=="College+"),ever_smoked=cfmt(d$smoking==1),alcohol_last_year=cfmt(d$alcohol==1),income_asinh=qfmt(d$income_asinh),asset_count=qfmt(d$asset_count),HD_z=qfmt(d$HD_z),KDM_BAA_z=qfmt(d$KDM_BAA_z))}))
write.csv(tab1,file.path(outdir,"table1.csv"),row.names=FALSE)

# Sensitivities: alternate HD/KDM, cluster level, primary-fuel-only trajectory,
# exclusion of raw biomarker extremes, and complete covariate records.
sens<-models %>% filter(analysis %in% c("Fully adjusted","Fully adjusted; household clustered","Overlap weighted"))
sens<-bind_rows(sens,std_contrast("HD_alt_z",full_rhs,primary,"Alternative HD: full-adult reference + 1/99 winsor"),std_contrast("KDM_BAA_direct_z",full_rhs,primary,"Alternative KDM: direct BA-age"))
extreme<-rep(FALSE,nrow(dat));for(j in names(Xraw)){qq<-quantile(Xraw[[j]][adult_complete],c(.005,.995),na.rm=TRUE);extreme<-extreme|(Xraw[[j]]<qq[1]|Xraw[[j]]>qq[2])}
analysis$any_extreme<-extreme[eligible]; p_noext<-analysis%>%filter(trajectory%in%c("persistent_any_solid","solid_to_clean_only"),!any_extreme)%>%droplevels()
sens<-bind_rows(sens,std_contrast("HD_z",full_rhs,p_noext,"Exclude any biomarker outside 0.5/99.5%"),std_contrast("KDM_BAA_z",full_rhs,p_noext,"Exclude any biomarker outside 0.5/99.5%"))
ccvars<-c("education_level_code","hhincpc_cpi","asset_count","smoking","alcohol");pcc<-primary[complete.cases(primary[,ccvars]),]
sens<-bind_rows(sens,std_contrast("HD_z",full_rhs,pcc,"Complete covariates only"),std_contrast("KDM_BAA_z",full_rhs,pcc,"Complete covariates only"))

alt <- dat %>% filter(age>=18,!is.na(sex),known_primary_2009,n_prior_primary>=2,!is.na(HD_z),!is.na(KDM_BAA_z),trajectory_primary_only%in%c("persistent_any_solid","solid_to_clean_only"))
alt$trajectory<-factor(alt$trajectory_primary_only,levels=c("persistent_any_solid","solid_to_clean_only"))
sens<-bind_rows(sens,std_contrast("HD_z",full_rhs,alt,"Primary fuel only (stacking ignored)"),std_contrast("KDM_BAA_z",full_rhs,alt,"Primary fuel only (stacking ignored)"))
write.csv(sens,file.path(outdir,"sensitivity_results.csv"),row.names=FALSE)

# Missing-data feasibility (MI is not used automatically because missing-category
# coding preserves the target sample; report whether a future MAR analysis is viable).
mi_feas <- miss %>% mutate(mi_candidate = variable %in% c("education","income_asinh","asset_count","smoking_f","alcohol_f"), feasible = pct_missing < 30)
write.csv(mi_feas,file.path(outdir,"mi_feasibility.csv"),row.names=FALSE)

traj_counts<-analysis%>%count(trajectory,name="n");write.csv(traj_counts,file.path(outdir,"trajectory_counts_formal.csv"),row.names=FALSE)
outcome_contract<-data.frame(outcome=c("HD_z","KDM_BAA_z"),role=c("Primary","Secondary"),time="2009 cross-sectional",
 method=c("Log Mahalanobis distance; 20-39 y reference; locked 9-system biomarkers; standardized","5-fold cross-fitted sex-stratified KDM; fold-specific calibration and BAA residualization; standardized"),
 interpretation=c("Higher = greater multivariate physiological dysregulation","Higher = older biological profile relative to chronological age"))
write.csv(outcome_contract,file.path(outdir,"outcome_construction_contract.csv"),row.names=FALSE)

cat("Formal n:",nrow(analysis)," primary contrast n:",nrow(primary),"\n")
print(models)
cat("C2 formal analysis completed:",format(Sys.time()),"\n")
