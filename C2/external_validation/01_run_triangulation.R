options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(svglite); library(ragg)
})

set.seed(20260809)
wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(wd) != "external_validation") stop("Run from C2/external_validation")
ext_root <- Sys.getenv("C2_TRI_EXTERNAL_ROOT")
if (!nzchar(ext_root) || !dir.exists(ext_root)) stop("Set C2_TRI_EXTERNAL_ROOT")
dir.create("results", FALSE, TRUE); dir.create("figures", FALSE, TRUE); dir.create("logs", FALSE, TRUE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
num <- function(x) { x <- as.numeric(x); x[x < 0 | x == 999] <- NA_real_; x }
qfun <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

markers <- c("crp", "glu", "tg", "hdl", "wbc", "hgb", "sbp")
pretty_marker <- c(crp="hsCRP", glu="Glucose", tg="Triglycerides", hdl="HDL-C",
                   wbc="White blood cells", hgb="Haemoglobin", sbp="Mean SBP")
transform_marker <- function(d) d %>% mutate(
  crp=log1p(crp), glu=log(glu), tg=log(tg), hdl=hdl,
  wbc=log(wbc), hgb=hgb, sbp=sbp)

fit_frozen <- function(raw, version) {
  d0 <- raw %>% filter(age >= 18, gender %in% c(1,2))
  d <- transform_marker(d0)
  lim <- lapply(markers, function(v) c(lo=qfun(d[[v]], .005), hi=qfun(d[[v]], .995)))
  names(lim) <- markers
  for (v in markers) d[[v]] <- pmin(pmax(d[[v]], lim[[v]][1]), lim[[v]][2])
  cc <- complete.cases(d[, c("age","gender",markers)])
  d <- d[cc, ]
  young <- d %>% filter(age >= 20, age <= 39)
  if (nrow(young) < 200) stop("Insufficient CHNS young reference")
  center <- sapply(young[markers], median)
  cv <- cov(young[markers]); ridge <- max(diag(cv), na.rm=TRUE) * 1e-6
  inv <- solve(cv + diag(ridge, length(markers)))
  md2 <- rowSums(((as.matrix(d[markers]) - rep(center, each=nrow(d))) %*% inv) *
                   (as.matrix(d[markers]) - rep(center, each=nrow(d))))
  hd_raw <- log1p(sqrt(pmax(md2, 0)))
  hd_mean <- mean(hd_raw); hd_sd <- sd(hd_raw)
  sexpars <- list()
  for (sx in c(1,2)) {
    z <- d[d$gender == sx, ]
    fits <- lapply(markers, function(v) lm(z[[v]] ~ z$age))
    q <- sapply(fits, function(f) unname(coef(f)[1])); k <- sapply(fits, function(f) unname(coef(f)[2]))
    se <- sapply(fits, function(f) summary(f)$sigma)
    den <- sum(k^2/se^2)
    bae <- as.numeric((as.matrix(z[markers]) - matrix(q, nrow(z), length(q), byrow=TRUE)) %*% (k/se^2) / den)
    sba2 <- var(bae - z$age)
    ba <- (bae * den + z$age/sba2)/(den + 1/sba2)
    acc <- lm(ba ~ z$age)
    sexpars[[as.character(sx)]] <- list(q=q, k=k, se=se, den=den, sba2=sba2,
                                        acc=unname(coef(acc)), n=nrow(z))
  }
  list(version=version, markers=markers, limits=lim, center=center, inv=inv,
       hd_mean=hd_mean, hd_sd=hd_sd, sexpars=sexpars, n=nrow(d), young_n=nrow(young))
}

apply_frozen <- function(raw, par) {
  d <- transform_marker(raw)
  outside <- lapply(markers, function(v) {
    c(low=mean(d[[v]] < par$limits[[v]][1], na.rm=TRUE),
      high=mean(d[[v]] > par$limits[[v]][2], na.rm=TRUE))
  }); names(outside) <- markers
  for (v in markers) d[[v]] <- pmin(pmax(d[[v]], par$limits[[v]][1]), par$limits[[v]][2])
  d$complete <- complete.cases(d[, c("age","gender",markers)]) & d$gender %in% c(1,2)
  d$hd7_z <- d$kdm_ba <- d$kdm_baa <- NA_real_
  ii <- which(d$complete)
  X <- as.matrix(d[ii, markers]); delta <- X - rep(par$center, each=length(ii))
  md2 <- rowSums((delta %*% par$inv) * delta)
  d$hd7_z[ii] <- (log1p(sqrt(pmax(md2,0))) - par$hd_mean)/par$hd_sd
  for (sx in c(1,2)) {
    jj <- ii[d$gender[ii] == sx]; p <- par$sexpars[[as.character(sx)]]
    Xs <- as.matrix(d[jj, markers])
    bae <- as.numeric((Xs - matrix(p$q, nrow(Xs), length(p$q), byrow=TRUE)) %*% (p$k/p$se^2) / p$den)
    ba <- (bae*p$den + d$age[jj]/p$sba2)/(p$den + 1/p$sba2)
    d$kdm_ba[jj] <- ba
    d$kdm_baa[jj] <- ba - (p$acc[1] + p$acc[2]*d$age[jj])
  }
  attr(d, "outside") <- outside
  d
}

perf <- function(d) {
  z <- d %>% filter(complete, is.finite(kdm_ba), is.finite(age))
  fit <- lm(kdm_ba ~ age, z)
  data.frame(n=nrow(z), MAE=mean(abs(z$kdm_ba-z$age)), RMSE=sqrt(mean((z$kdm_ba-z$age)^2)),
             R2=cor(z$kdm_ba,z$age)^2, calibration_intercept=coef(fit)[1],
             calibration_slope=coef(fit)[2], BA_age_correlation=cor(z$kdm_ba,z$age),
             BAA_age_correlation=cor(z$kdm_baa,z$age), HD7_mean=mean(z$hd7_z), HD7_SD=sd(z$hd7_z))
}

# CHNS discovery source and frozen shared-seven-marker outcome models
chns0 <- readRDS("../data/selection_cohort_internal_restricted.rds")
chns_base <- chns0 %>% transmute(age=as.numeric(age), gender=as.numeric(gender),
  crp=as.numeric(HS_CRP), glu=as.numeric(glucose), tg=as.numeric(tg), hdl=as.numeric(HDL_C),
  wbc=as.numeric(wbc), hgb=as.numeric(hgb), sbp=as.numeric(sbp))
chns_hba <- chns0 %>% transmute(age=as.numeric(age), gender=as.numeric(gender),
  crp=as.numeric(HS_CRP), glu=as.numeric(HbA1c), tg=as.numeric(tg), hdl=as.numeric(HDL_C),
  wbc=as.numeric(wbc), hgb=as.numeric(hgb), sbp=as.numeric(sbp))
par_glu <- fit_frozen(chns_base, "glucose_shared7_primary")
par_hba <- fit_frozen(chns_hba, "hba1c_shared7_sensitivity")
saveRDS(list(glucose=par_glu, hba1c=par_hba), "results/chns_frozen_shared7_parameters.rds")
write_csv(data.frame(version=c(par_glu$version,par_hba$version), chns_complete_n=c(par_glu$n,par_hba$n),
                     young_reference_n=c(par_glu$young_n,par_hba$young_n)), "results/chns_algorithm_training_counts.csv")

# ELSA algorithm transport
elsa_dir <- file.path(ext_root, "ELSA_UK", "UKDA-5050-stata", "stata", "stata13_se")
gh <- read_dta(file.path(elsa_dir, "gh_elsa_h.dta"), col_select=c("idauniq","r8agey","ragender"))
nu <- read_dta(file.path(elsa_dir, "wave_8_elsa_nurse_data_eul_v1.dta"),
               col_select=c("idauniq","sys2","sys3","hdl","trig","hscrp","fglu","hba1c","hgb","wbc"))
elsa_join <- inner_join(gh, nu, by="idauniq")
elsa_raw <- elsa_join %>% transmute(age=num(r8agey), gender=num(ragender),
  crp=num(hscrp), glu=num(fglu), hba_ifcc=num(hba1c), tg=num(trig), hdl=num(hdl),
  wbc=num(wbc), hgb=num(hgb), s2=num(sys2),s3=num(sys3)) %>%
  mutate(sbp=ifelse(is.finite(s2)&is.finite(s3),(s2+s3)/2,NA_real_)) %>% select(-s2,-s3)
elsa_hba <- elsa_raw %>% mutate(glu=0.09148*hba_ifcc + 2.152)
elsa_g <- apply_frozen(select(elsa_raw,-hba_ifcc), par_glu)
elsa_h <- apply_frozen(select(elsa_hba,-hba_ifcc), par_hba)
pg <- perf(elsa_g) %>% mutate(version="glucose_shared7_primary", .before=1)
ph <- perf(elsa_h) %>% mutate(version="hba1c_shared7_sensitivity", .before=1)
transport <- bind_rows(pg,ph) %>% mutate(
  pure_transport_pass=MAE<=15 & R2>=.40 & calibration_slope>=.50 & calibration_slope<=1.50 & abs(BAA_age_correlation)<.10)
primary_pass <- transport$pure_transport_pass[transport$version=="glucose_shared7_primary"]
write_csv(transport, "results/elsa_transport_performance.csv")

flow <- data.frame(stage=c("Gateway wave-8 records","Nurse wave-8 records","Linked records",
                           "Complete glucose shared-seven","Complete HbA1c shared-seven"),
                   n=c(nrow(gh),nrow(nu),nrow(elsa_join),sum(elsa_g$complete),sum(elsa_h$complete)))
write_csv(flow, "results/elsa_sample_flow.csv")
miss <- data.frame(marker=c("Age","Sex",pretty_marker,"HbA1c (IFCC)"),
  available_n=c(sum(!is.na(elsa_raw$age)),sum(!is.na(elsa_raw$gender)),sapply(elsa_raw[markers],function(x)sum(!is.na(x))),sum(!is.na(elsa_raw$hba_ifcc))))
miss$missing_n <- nrow(elsa_raw)-miss$available_n; miss$missing_pct <- 100*miss$missing_n/nrow(elsa_raw)
write_csv(miss, "results/elsa_marker_missingness.csv")

shift_one <- function(eraw, craw, par, ver) bind_rows(lapply(markers,function(v){
  et <- transform_marker(eraw)[[v]]; ct <- transform_marker(craw)[[v]]
  data.frame(version=ver,marker=pretty_marker[[v]],elsa_n=sum(is.finite(et)),
    chns_median=median(ct,na.rm=TRUE),elsa_median=median(et,na.rm=TRUE),
    below_CHNS_winsor_pct=100*mean(et<par$limits[[v]][1],na.rm=TRUE),
    above_CHNS_winsor_pct=100*mean(et>par$limits[[v]][2],na.rm=TRUE))
}))
write_csv(bind_rows(shift_one(select(elsa_raw,-hba_ifcc),chns_base,par_glu,par_glu$version),
                    shift_one(select(elsa_hba,-hba_ifcc),chns_hba,par_hba,par_hba$version)),
          "results/elsa_marker_shift_diagnostics.csv")

boot_metrics <- function(d, ver, B=1000) {
  z <- d %>% filter(complete)
  bind_rows(lapply(seq_len(B),function(b){
    x <- z[sample.int(nrow(z),nrow(z),replace=TRUE),]
    cbind(data.frame(version=ver,replicate=b),perf(x))
  }))
}
bg <- boot_metrics(elsa_g,par_glu$version); bh <- boot_metrics(elsa_h,par_hba$version)
boots <- bind_rows(bg,bh); write_csv(boots,"results/elsa_bootstrap_replicates.csv")
metric_names <- c("MAE","RMSE","R2","calibration_intercept","calibration_slope","BAA_age_correlation","HD7_mean")
boot_long <- boots %>% pivot_longer(all_of(metric_names),names_to="metric",values_to="estimate")
boot_sum <- boot_long %>% group_by(version,metric) %>% summarise(
  lower=qfun(estimate,.025),upper=qfun(estimate,.975),estimate=median(estimate),.groups="drop")
write_csv(boot_sum,"results/elsa_bootstrap_summary.csv")
elsa_fig <- elsa_g %>% filter(complete) %>% transmute(study_id=row_number(),age,kdm_ba,kdm_baa,hd7_z,gender)
write_csv(elsa_fig,"results/elsa_primary_deidentified.csv")

# MHAS six-wave trajectory-structure replication
mhas <- read_dta(file.path(ext_root,"MHAS_Mexico","H_MHAS_d.dta"),
                 col_select=c("unhhidnp",paste0("hh",1:6,"cookfuel_m")))
fuel <- mhas %>% transmute(across(all_of(paste0("hh",1:6,"cookfuel_m")),~case_when(as.numeric(.x)==1~"gas",as.numeric(.x)==2~"solid",TRUE~NA_character_)))
names(fuel) <- paste0("w",1:6)
pre_n <- rowSums(!is.na(fuel[paste0("w",1:5)])); eligible <- !is.na(fuel$w6) & pre_n>=2
classify <- function(x) {
  y <- x[!is.na(x)]; if(length(y)<3) return("ineligible")
  if(all(y=="gas")) return("persistent gas")
  if(all(y=="solid")) return("persistent solid")
  code <- ifelse(y=="gas",1,0); dif <- diff(code)
  if(y[1]=="solid" && tail(y,1)=="gas" && all(dif>=0)) return("solid-to-gas")
  if(y[1]=="gas" && tail(y,1)=="solid") return("gas-to-solid")
  "complex"
}
traj <- rep("ineligible",nrow(fuel)); traj[eligible] <- apply(fuel[eligible,],1,classify)
lev <- c("persistent solid","solid-to-gas","persistent gas","gas-to-solid","complex","ineligible")
tc <- as.data.frame(table(factor(traj,levels=lev)),responseName="n")
names(tc)[1] <- "trajectory"; tc$eligible_denominator <- sum(eligible)
tc$percent_among_eligible <- ifelse(tc$trajectory=="ineligible",NA,100*tc$n/sum(eligible))
write_csv(tc,"results/mhas_trajectory_counts.csv")
wave_dist <- bind_rows(lapply(1:6,function(w){
  x <- fuel[[w]]; data.frame(wave=w,state=c("gas","solid","unknown"),
    n=c(sum(x=="gas",na.rm=TRUE),sum(x=="solid",na.rm=TRUE),sum(is.na(x))))
})) %>% group_by(wave) %>% mutate(percent=100*n/sum(n)) %>% ungroup()
write_csv(wave_dist,"results/mhas_wave_distribution.csv")
trans <- bind_rows(lapply(1:5,function(w){
  x <- fuel[[w]]; y <- fuel[[w+1]]; ok=!is.na(x)&!is.na(y)
  as.data.frame(table(from=factor(x[ok],levels=c("solid","gas")),to=factor(y[ok],levels=c("solid","gas")))) %>%
    mutate(from_wave=w,to_wave=w+1)
}))
write_csv(trans,"results/mhas_transition_matrix.csv")
hist <- as.data.frame(table(valid_pre_waves=pre_n),responseName="n") %>% mutate(percent=100*n/sum(n))
write_csv(hist,"results/mhas_history_length.csv")

# Nature-style R-only figure
theme_nat <- theme_classic(base_family="Arial",base_size=8) + theme(
  plot.title=element_text(face="bold",size=9,hjust=0),plot.subtitle=element_text(size=7,color="#444444"),
  axis.title=element_text(size=7.5),axis.text=element_text(size=6.8,color="#222222"),
  strip.background=element_blank(),strip.text=element_text(face="bold",size=7),
  plot.margin=margin(5,7,5,5),legend.position="none")
pal <- c(CHNS="#0072B2",ELSA="#009E73",MHAS="#E69F00")

pa_dat <- data.frame(x=c(0,1,-1),y=c(1,-.48,-.48),cohort=c("CHNS","ELSA","MHAS"),
                     role=c("Discovery association","Outcome-algorithm transport","Fuel-trajectory structure"),
                     role_y=c(.48,-1.02,-1.02))
pa <- ggplot(pa_dat,aes(x,y))+
  annotate("segment",x=0,y=.78,xend=.79,yend=-.32,arrow=arrow(length=unit(2,"mm")),color="#777777",linewidth=.45)+
  annotate("segment",x=0,y=.78,xend=-.79,yend=-.32,arrow=arrow(length=unit(2,"mm")),color="#777777",linewidth=.45)+
  geom_point(aes(fill=cohort),shape=21,size=12,color="white",stroke=.7)+scale_fill_manual(values=pal)+
  geom_text(aes(label=cohort),color="white",fontface="bold",size=2.7)+
  geom_text(aes(y=role_y,label=role),size=2.25,lineheight=.92)+
  annotate("text",x=0,y=-1.42,label="Independent modules — no pooled effect",size=2.25,color="#555555")+
  coord_cartesian(xlim=c(-1.5,1.5),ylim=c(-1.58,1.28),clip="off")+theme_void(base_family="Arial")+
  theme(legend.position="none")+
  ggtitle("a  Triangular validation design")+theme(plot.title=element_text(face="bold",size=9,hjust=0),plot.margin=margin(5,7,5,5))

cal <- coef(lm(kdm_ba~age,elsa_fig))
pb <- ggplot(elsa_fig,aes(age,kdm_ba))+geom_point(alpha=.18,size=.65,color=pal[["ELSA"]])+
  geom_abline(slope=1,intercept=0,linetype=2,color="#777777",linewidth=.45)+
  geom_abline(slope=cal[2],intercept=cal[1],color="#D55E00",linewidth=.65)+
  annotate("text",x=Inf,y=-Inf,hjust=1.02,vjust=-.5,size=2.25,
           label=sprintf("n=%s; slope=%.2f",format(nrow(elsa_fig),big.mark=","),cal[2]))+
  labs(title="b  ELSA frozen-algorithm calibration",subtitle="CHNS shared-seven-marker KDM applied without refitting",x="Chronological age (years)",y="Predicted biological age (years)")+theme_nat

show_metrics <- c(MAE="MAE (years)",R2="R²",calibration_slope="Calibration slope",BAA_age_correlation="BAA–age correlation")
pcd <- boot_sum %>% filter(version==par_glu$version,metric%in%names(show_metrics)) %>%
  mutate(label=unname(show_metrics[metric]),target=case_when(metric=="MAE"~15,metric=="R2"~.40,metric=="calibration_slope"~1,TRUE~0))
pc <- ggplot(pcd,aes(estimate,label))+geom_vline(aes(xintercept=target),linetype=2,color="#888888",linewidth=.4)+
  geom_errorbar(aes(xmin=lower,xmax=upper),width=.16,orientation="y",color=pal[["ELSA"]],linewidth=.55)+
  geom_point(size=2.2,shape=21,fill=pal[["ELSA"]],color="white",stroke=.35)+
  facet_wrap(~label,scales="free",ncol=2)+labs(title="c  ELSA transport diagnostics",subtitle="Points and 95% participant-bootstrap intervals",x=NULL,y=NULL)+
  theme_nat+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

td <- tc %>% filter(trajectory!="ineligible") %>% mutate(trajectory=factor(trajectory,levels=rev(lev[1:5])))
pd <- ggplot(td,aes(n,trajectory,fill=trajectory))+geom_col(width=.62)+
  geom_text(aes(label=sprintf("%s (%.1f%%)",format(n,big.mark=","),percent_among_eligible)),hjust=-.05,size=2.3)+
  scale_fill_manual(values=c("persistent solid"="#D55E00","solid-to-gas"="#E69F00","persistent gas"="#009E73","gas-to-solid"="#56B4E9","complex"="#999999"))+
  scale_x_continuous(expand=expansion(mult=c(0,.24)))+
  labs(title="d  MHAS fuel-trajectory replication",subtitle=sprintf("Six waves; eligible n=%s",format(sum(eligible),big.mark=",")),x="Participants",y=NULL)+theme_nat

fig <- (pa|pb)/(pc|pd)+plot_layout(heights=c(1,1.03))
ggsave("figures/Figure3_C2_triangular_validation.svg",fig,width=183,height=130,units="mm",device=svglite)
ggsave("figures/Figure3_C2_triangular_validation.pdf",fig,width=183,height=130,units="mm",device=cairo_pdf)
agg_tiff("figures/Figure3_C2_triangular_validation.tiff",width=4323,height=3071,units="px",res=600,compression="lzw"); print(fig); dev.off()
agg_png("figures/Figure3_C2_triangular_validation.png",width=2161,height=1535,units="px",res=300); print(fig); dev.off()
write_csv(pa_dat,"figures/Figure3_source_panel_a.csv")
write_csv(elsa_fig,"figures/Figure3_source_panel_b.csv")
write_csv(pcd,"figures/Figure3_source_panel_c.csv")
write_csv(td,"figures/Figure3_source_panel_d.csv")

cat("TRIANGULATION_RUN=PASS\n")
cat(sprintf("ELSA_PRIMARY_N=%d PURE_TRANSPORT_PASS=%s\n",sum(elsa_g$complete),primary_pass))
cat(sprintf("MHAS_ELIGIBLE_N=%d STRUCTURE_PASS=%s\n",sum(eligible),all(tc$n[match(c("persistent solid","solid-to-gas"),tc$trajectory)]>=200)))
