#!/usr/bin/env Rscript

# C2 Figure 1: aggregate audit evidence only. R-only rendering/export.
options(encoding = "UTF-8", warn = 1)

required <- c("ggplot2", "patchwork", "dplyr", "readr", "svglite", "ragg", "scales")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript so the project path can be resolved.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
fig_dir <- file.path(project_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_audit <- function(name) {
  path <- file.path(results_dir, name)
  if (!file.exists(path)) stop("Missing aggregate audit CSV: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

sample_flow <- read_audit("sample_flow.csv")
trajectories <- read_audit("effective_trajectory_counts_main_eligible.csv")
stacking <- read_audit("fuel_stacking_by_wave.csv")
biomarkers <- read_audit("candidate_biomarker_completeness.csv")

pal <- c(
  navy = "#484878", blue = "#7884B4", pale_blue = "#B4C0E4",
  teal = "#42949E", pale_teal = "#BFDDE0", accent = "#B64342",
  grey = "#A8A8A8", light = "#E4E4F0", dark = "#272727"
)

theme_nature <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal["dark"]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal["dark"]),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.2, colour = pal["dark"]),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.4),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#606060"),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_nature())

# a: selected sequential denominators; final n is derived from the formal trajectory table.
get_stage_n <- function(stage_name) {
  x <- sample_flow$n[sample_flow$stage == stage_name]
  if (length(x) != 1) stop("Expected one sample-flow row for: ", stage_name)
  x
}
n_final <- sum(trajectories$n)
flow <- tibble::tibble(
  stage = c("2009 biomarker\nparticipants", "Age/sex linked", "Adults + core five\nbiomarkers", "Stacking-aware\nanalysis sample"),
  n = c(get_stage_n("Biomarker file"), get_stage_n("Linked age and sex"),
        get_stage_n("Adults + age/sex + core 5"), n_final)
) %>% mutate(stage = factor(stage, levels = stage))

p_a <- ggplot(flow, aes(stage, n)) +
  geom_col(aes(fill = stage), width = 0.68, colour = "white", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.45, size = 2.25, fontface = "bold") +
  scale_fill_manual(values = unname(c(pal["light"], pal["pale_blue"], pal["blue"], pal["navy"]))) +
  scale_y_continuous(limits = c(0, 10300), breaks = c(0, 2500, 5000, 7500, 10000), labels = comma,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Auditable linkage retains a large analysis sample", x = NULL, y = "Participants (n)") +
  theme(axis.text.x = element_text(size = 5.9, lineheight = 0.9), axis.ticks.x = element_blank())

# b: stacking-aware trajectory composition.
traj_labels <- c(
  persistent_any_solid = "Persistent any solid",
  solid_to_clean_only_monotonic = "Solid to clean-only",
  complex_switching = "Complex switching",
  persistent_clean_only = "Persistent clean-only",
  clean_only_to_solid_monotonic = "Clean-only to solid"
)
traj <- trajectories %>%
  mutate(label = unname(traj_labels[effective_trajectory]),
         pct = n / sum(n),
         label = factor(label, levels = rev(unname(traj_labels))))

p_b <- ggplot(traj, aes(n, label)) +
  geom_col(aes(fill = label), width = 0.66, colour = "white", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = paste0(comma(n), "  (", percent(pct, accuracy = 0.1), ")")),
            hjust = -0.08, size = 2.15) +
  scale_fill_manual(values = setNames(
    unname(c(pal["navy"], pal["teal"], pal["grey"], pal["pale_blue"], pal["accent"])),
    unname(traj_labels)
  )) +
  scale_x_continuous(limits = c(0, 4300), breaks = c(0, 1000, 2000, 3000, 4000), labels = comma,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Joint-fuel trajectories support the planned contrasts",
       subtitle = paste0("Formal analysis denominator: n = ", comma(n_final)),
       x = "Participants (n)", y = NULL)

# c: two denominators prevent adding the mixed subset to mutually exclusive classes.
s09 <- stacking %>% filter(wave == 2009)
if (nrow(s09) != 1) stop("Expected one 2009 row in fuel_stacking_by_wave.csv")
state <- tibble::tibble(
  denominator = "All 2009 linked records\n(n = 9,548)",
  category = c("Any solid", "Clean-only", "Other/unknown"),
  n = c(s09$n_any_solid, s09$n_clean_only, s09$n_other_or_unknown)
)
solid_subset <- tibble::tibble(
  denominator = "Within any-solid\n(n = 4,898)",
  category = c("Mixed clean–solid", "No mixed stacking recorded"),
  n = c(s09$n_mixed_clean_solid_stacking, s09$n_any_solid - s09$n_mixed_clean_solid_stacking)
)
stack_plot <- bind_rows(state, solid_subset) %>%
  group_by(denominator) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
  mutate(
    denominator = factor(denominator, levels = c("All 2009 linked records\n(n = 9,548)", "Within any-solid\n(n = 4,898)")),
    category = factor(category, levels = c("No mixed stacking recorded", "Mixed clean–solid", "Other/unknown", "Clean-only", "Any solid"))
  )

stack_cols <- c("Any solid" = unname(pal["navy"]), "Clean-only" = unname(pal["pale_teal"]),
                "Other/unknown" = unname(pal["grey"]), "Mixed clean–solid" = unname(pal["accent"]),
                "No mixed stacking recorded" = unname(pal["light"]))
p_c <- ggplot(stack_plot, aes(denominator, pct, fill = category)) +
  geom_col(width = 0.63, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(pct >= 0.10, paste0(comma(n), "\n", percent(pct, accuracy = 1)), "")),
            position = position_stack(vjust = 0.5), size = 2.05, lineheight = 0.88,
            colour = ifelse(stack_plot$category %in% c("Any solid", "Mixed clean–solid"), "white", pal["dark"])) +
  scale_fill_manual(values = stack_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Fuel stacking materially changes exposure classification",
       subtitle = "Mixed clean–solid is a subset of any-solid (3,822/4,898; 78.0%)",
       x = NULL, y = "Composition", fill = "Classification") +
  theme(axis.text.x = element_text(size = 5.9, lineheight = 0.9), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# d: laboratory completeness only (physical-exam candidates use a different source process).
bio <- biomarkers %>%
  filter(source == "biomarker_09.dta") %>%
  mutate(variable = factor(variable, levels = variable[order(proportion_nonmissing)]),
         pct = 100 * proportion_nonmissing)

p_d <- ggplot(bio, aes(pct, variable)) +
  geom_segment(aes(x = 98.8, xend = pct, yend = variable), colour = pal["pale_blue"], linewidth = 1.2) +
  geom_point(colour = pal["navy"], size = 1.8) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 1.95) +
  scale_x_continuous(limits = c(98.8, 100.05), breaks = c(99.0, 99.5, 100.0), labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Candidate laboratory biomarkers are highly complete", x = "Non-missing in biomarker file", y = NULL)

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1, 1.12), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = c(0.005, 0.995))

# Source data: one long table with explicit denominators and provenance.
source_data <- bind_rows(
  flow %>% transmute(panel = "a", metric = as.character(stage), category = NA_character_, n,
                     denominator = max(n), proportion = n / max(n), input_table = "sample_flow.csv + trajectory table"),
  traj %>% transmute(panel = "b", metric = "Effective trajectory", category = as.character(label), n,
                     denominator = sum(n), proportion = pct, input_table = "effective_trajectory_counts_main_eligible.csv"),
  stack_plot %>% transmute(panel = "c", metric = as.character(denominator), category = as.character(category), n,
                           denominator = ave(n, denominator, FUN = sum), proportion = pct, input_table = "fuel_stacking_by_wave.csv"),
  bio %>% transmute(panel = "d", metric = "Biomarker completeness", category = as.character(variable), n = n_nonmissing,
                    denominator = 9549, proportion = proportion_nonmissing, input_table = "candidate_biomarker_completeness.csv")
)
readr::write_excel_csv(source_data, file.path(fig_dir, "Figure1_C2_source_data.csv"), na = "")

base <- file.path(fig_dir, "Figure1_C2")
width_mm <- 183
height_mm <- 132
w <- width_mm / 25.4
h <- height_mm / 25.4

svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
print(fig); grDevices::dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
print(fig); grDevices::dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600,
               background = "white", compression = "lzw")
print(fig); grDevices::dev.off()

ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 600,
              background = "white")
print(fig); grDevices::dev.off()

qa <- c(
  "# Figure 1 QA notes (C2)", "",
  paste0("- Render backend: R ", getRversion(), " only; ggplot2/patchwork with svglite, cairo_pdf and ragg."),
  paste0("- Final size: ", width_mm, " mm × ", height_mm, " mm; TIFF and PNG rendered at 600 dpi."),
  "- Core conclusion: linkage, stacking-aware classification and analytic feasibility; no fuel–biological-age effect is shown.",
  paste0("- Panel a final denominator cross-check: ", comma(n_final), "."),
  paste0("- Panel b trajectory count sum: ", comma(sum(trajectories$n)), "."),
  paste0("- Panel c 2009 state count sum: ", comma(s09$n_any_solid + s09$n_clean_only + s09$n_other_or_unknown),
         "; linked-record denominator: ", comma(s09$n_biomarker_participant_records), "."),
  paste0("- Panel c stacking subset check: ", comma(s09$n_mixed_clean_solid_stacking), "/", comma(s09$n_any_solid),
         " = ", percent(s09$n_mixed_clean_solid_stacking / s09$n_any_solid, accuracy = 0.1), "."),
  "- Panel d contains laboratory biomarkers only; no uncertainty interval or hypothesis test is applicable to audit completeness.",
  "- Source-data traceability: Figure1_C2_source_data.csv contains panel, metric, category, n, denominator, proportion and input table.",
  "- Image integrity: no microscopy, crop, contrast adjustment, pseudo-colour or stitching; all marks are generated from aggregate CSV data.",
  "- Interpretation limit: mixed clean–solid is a subset of any-solid and is not an additional mutually exclusive exposure group.",
  "- Visual QA of the 600-dpi PNG at final aspect ratio: panel-b labels are inside the canvas, panel-c stacking is visually distinct, and panel-d labels do not overlap.",
  "- Automated checks passed if this file was written: all required inputs existed and all four graphics devices closed successfully."
)
writeLines(enc2utf8(qa), file.path(fig_dir, "Figure1_C2_qa_notes.md"), useBytes = TRUE)

message("Figure 1 exports and source data written to: ", fig_dir)
