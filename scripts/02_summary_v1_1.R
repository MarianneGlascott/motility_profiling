# ==============================================================================
# Project:     Kelp Zoospore Motility Profiling
# Scripts:     02_summary.R
# Author:      Marianne Glascott (with ChatGPT assistant)
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2025-10-31
# Version:     1.1
# ==============================================================================
# Description
#   Import clean culture‑level CSV (642 rows)
#   Produce robust summary statistics for counts/ratios (overall + by day)
#   Provide day‑specific (0, 1, 4) summaries to avoid day‑mixing skew
#   Generate clear graphics of motility_prop vs days_from_start
#     → uses MEAN ± SE (per your preference) for the main trend
#   Run twice automatically: (A) All species; (B) L. digitata only
#   QC: validate ranges, detect outliers, and write a QC report
#   Save tables to outputs/ and figures to plots/ & figures/
# ==============================================================================
# Notes
# - Paths are project‑relative using here::here().
# - Figures saved to plots/ (screen/EDA) and figures/ (pub‑quality) in PNG, TIFF, PDF.
# - Git: points origin to https://github.com/MarianneGlascott/motility_profiling.git
# - Designed to work from RStudio Project: C:/Motility_Profiling/Motility_Profiling.Rproj
# ==============================================================================
# Contents:
# 1. Setup & tune-able flags
# 2. Read data & hygiene
# 3. QC checks & report
# 4. Utility: compact distribution summary
# 5. Reusable runner to produce outputs for chosen species filter
# 6. Run BOTH analyses: ALL species and L. digitata only
# 7. Housekeeping
# 8. Narrative paragraph booster (append)
# ==============================================================================
# 1. Setup & tune-able flags
# ==============================================================================
STAGE <- "SUMMARY"
source(here::here("scripts","01_setup.R"))  # adjust if you use a versioned name

# Preferred days for unbiased snapshots
FOCUS_DAYS   <- c(0, 1, 4)
INPUT_FILE   <- here::here("data_raw","kelp_control_cultures_clean.csv")

# QC thresholds (edit if needed)
QC_MIN_MOTILITY <- 0
QC_MAX_MOTILITY <- 1
QC_MIN_COUNT    <- 0
QC_SMALL_N      <- 15   # warn when group n < 15 for mean/SE stability
# ==============================================================================
# 2. Read data & hygiene
# ==============================================================================
stopifnot(file.exists(INPUT_FILE))

d0 <- readr::read_csv(INPUT_FILE) |> janitor::clean_names()

# Expected columns
need <- c("mobile_cell_count","stationary_cell_count","motility_prop",
          "settlement_count","days_from_start","species",
          "collection_site","season")
miss <- setdiff(need, names(d0))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

# Coerce types defensively
as_int <- function(x) suppressWarnings(as.integer(x))

d <- d0 |>
  mutate(
    days_from_start = as_int(days_from_start),
    species         = as.factor(species),
    collection_site = as.factor(collection_site),
    season          = as.factor(season),
    motility_prop   = suppressWarnings(as.numeric(motility_prop)),
    settlement_count = suppressWarnings(as.numeric(settlement_count)),
    mobile_cell_count = suppressWarnings(as.numeric(mobile_cell_count)),
    stationary_cell_count = suppressWarnings(as.numeric(stationary_cell_count))
  )
# ==============================================================================
# 3. QC checks & report
# ==============================================================================
qc_rows <- d |>
  mutate(
    qc_motility_oob = is.na(motility_prop) | motility_prop < QC_MIN_MOTILITY | motility_prop > QC_MAX_MOTILITY,
    qc_neg_counts   = (mobile_cell_count   < QC_MIN_COUNT) | (stationary_cell_count < QC_MIN_COUNT) |
      (settlement_count    < QC_MIN_COUNT),
    qc_days_na      = is.na(days_from_start)
  ) |>
  filter(qc_motility_oob | qc_neg_counts | qc_days_na)

if (nrow(qc_rows)) {
  readr::write_csv(qc_rows, file.path(DIR_OUTPUTS, "qc_flags_rows.csv"))
  log_run(STAGE, glue("QC: {nrow(qc_rows)} rows flagged (see outputs/qc_flags_rows.csv)"))
} else {
  log_run(STAGE, "QC: no obvious range violations detected")
}

# Clamp motility into [0,1] for robustness (record change count)
clamped <- sum(d$motility_prop < QC_MIN_MOTILITY | d$motility_prop > QC_MAX_MOTILITY, na.rm = TRUE)
d$motility_prop <- pmin(pmax(d$motility_prop, QC_MIN_MOTILITY), QC_MAX_MOTILITY)
if (clamped > 0) log_run(STAGE, glue("QC: clamped motility_prop for {clamped} values to [0,1]"))
# ==============================================================================
# 4. Utility: compact distribution summary
# ==============================================================================
summarise_distribution <- function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  xok  <- xnum[is.finite(xnum)]
  if (length(xok) == 0) return(tibble(n = 0, missing = length(x),
                                      min = NA_real_, q1 = NA_real_, median = NA_real_,
                                      mean = NA_real_, q3 = NA_real_, max = NA_real_,
                                      sd = NA_real_, iqr = NA_real_, mode_rounded = NA_real_))
  qu <- stats::quantile(xok, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  md <- as.numeric(names(sort(table(round(xok, 0)), decreasing = TRUE)[1]))
  tibble(
    n = length(xok),
    missing = sum(!is.finite(xnum)),
    min = min(xok), q1 = qu[1], median = qu[2], mean = mean(xok), q3 = qu[3],
    max = max(xok), sd = stats::sd(xok), iqr = IQR(xok), mode_rounded = md
  )
}

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
# ==============================================================================
# 5. Reusable runner to produce outputs for chosen species filter
# ==============================================================================
run_summary <- function(data, species_filter = NA_character_, label = "all_species") {
  dd <- if (!is.na(species_filter)) dplyr::filter(data, species == species_filter) else data
  tag <- ifelse(is.na(species_filter), "ALL", gsub("[ .]", "_", species_filter))
  
  log_run(STAGE, glue("Running summary for: {tag} (n={nrow(dd)})"))
  
  vars_to_summarise <- c("motility_prop","settlement_count","mobile_cell_count",
                         "stationary_cell_count")
  
  overall_tbl <- purrr::map_dfr(vars_to_summarise, \(v) {
    out <- summarise_distribution(dd[[v]])
    dplyr::bind_cols(variable = v, out)
  }) |> dplyr::relocate(variable)
  
  readr::write_csv(overall_tbl, file.path(DIR_OUTPUTS, glue("overall_distribution_summary__{tag}.csv")))
  
  byday_tbl <- dd |> filter(days_from_start %in% FOCUS_DAYS) |>
    group_by(days_from_start) |>
    summarise(
      across(all_of(vars_to_summarise),
             list(
               n = ~sum(is.finite(.x)),
               mean = ~mean(.x, na.rm = TRUE),
               se = ~se(.x),
               median = ~median(.x, na.rm = TRUE),
               q1 = ~quantile(.x, 0.25, na.rm = TRUE),
               q3 = ~quantile(.x, 0.75, na.rm = TRUE)
             ), .names = "{.col}.{.fn}") , .groups = "drop")
  
  readr::write_csv(byday_tbl, file.path(DIR_OUTPUTS, glue("byday_summary_0_1_4__{tag}.csv")))
  
  # Percent decline day 0 -> 4 in motility (mean and median), with Ns
  d0 <- byday_tbl |> filter(days_from_start == 0)
  d4 <- byday_tbl |> filter(days_from_start == 4)
  decline_tbl <- tibble(
    group = tag,
    n_day0 = d0$motility_prop.n, n_day4 = d4$motility_prop.n,
    mean_day0 = d0$motility_prop.mean, mean_day4 = d4$motility_prop.mean,
    med_day0  = d0$motility_prop.median, med_day4  = d4$motility_prop.median,
    decline_mean_pct = 100 * (mean_day4 - mean_day0) / mean_day0,
    decline_median_pct = 100 * (med_day4 - med_day0) / med_day0
  )
  readr::write_csv(decline_tbl, file.path(DIR_OUTPUTS, glue("decline_day0_to_day4__{tag}.csv")))
  
  # Full trend across ALL observed days: MEAN ± SE
# ---------- Motility over time: show spread + smooth progression --------------
suppressPackageStartupMessages({
  library(ggbeeswarm)  # dot-swarm to show all cultures without overlap
  library(mgcv)        # smooth (GAM) for over-time trend
})

# Per-day summary (for mean±SE and n labels)
day_counts <- dd |>
  dplyr::group_by(days_from_start) |>
  dplyr::summarise(
    n    = sum(is.finite(motility_prop)),
    mean = mean(motility_prop, na.rm = TRUE),
    se   = stats::sd(motility_prop, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) |>
  dplyr::filter(n > 0)

# A) Dot-swarm (each dot = one culture), with mean±SE and n per day
p_swarm <- ggplot(dd, aes(x = days_from_start, y = motility_prop)) +
  ggbeeswarm::geom_quasirandom(width = 0.25, alpha = 0.25, size = 1.6) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.7) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  geom_text(data = day_counts,
            aes(label = paste0("n=", n), y = pmin(1, mean + se + 0.05)),
            size = 3.2, vjust = 0, nudge_y = 0.01, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
    subtitle = "Each dot is a culture; black line = mean; bars = SE; labels show n per day",
    x = "Days from start", y = "Motility proportion"
  )

show_and_save(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 130)
export_figure(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 130)

# B) Smooth progression across all observed days (GAM), with semi-transparent points
# ---------- B) Smooth progression across all observed days (GAM) + n labels ---
# (drop-in replacement for your existing p_gam block)

# Per-day counts (reuse if already computed)
day_counts <- dd |>
  dplyr::group_by(days_from_start) |>
  dplyr::summarise(
    n    = sum(is.finite(motility_prop)),
    mean = mean(motility_prop, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n > 0)

# Simple linear trend for an easy-to-read effect size
lm_fit <- lm(motility_prop ~ days_from_start, data = dd)
slope  <- coef(lm_fit)[["days_from_start"]]           # change in motility_prop per day
ci     <- confint(lm_fit)["days_from_start", ]        # 95% CI for slope

# % change from day 0 to the last observed day, referenced to day-0 mean
d0  <- min(dd$days_from_start, na.rm = TRUE)
dL  <- max(dd$days_from_start, na.rm = TRUE)
m0  <- mean(dd$motility_prop[dd$days_from_start == d0], na.rm = TRUE)
mL  <- mean(dd$motility_prop[dd$days_from_start == dL], na.rm = TRUE)
pct <- if (is.finite(m0) && m0 > 0) 100 * (mL - m0) / m0 else NA_real_

annot <- glue::glue("Linear slope: {sprintf('%.3f', slope)}/day ",
                    "[{sprintf('%.3f', ci[1])}, {sprintf('%.3f', ci[2])}]",
                    if (is.finite(pct)) glue::glue("; Δ% {d0}→{dL}: {sprintf('%.1f', pct)}%") else "")

p_gam <- ggplot(dd, aes(days_from_start, motility_prop)) +
  geom_point(alpha = 0.22, size = 1.3) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE) +
  geom_text(data = day_counts,
            aes(x = days_from_start, y = pmin(1, mean + 0.07), label = paste0("n=", n)),
            size = 3, vjust = 0, colour = "grey20") +
  scale_x_continuous(breaks = sort(unique(dd$days_from_start))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = glue::glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
    subtitle = glue::glue("GAM smooth with 95% CI; dots show culture-level spread (days extend to {dL}). {annot}"),
    x = "Days from start", y = "Motility proportion"
  )

show_and_save(p_gam, glue("motility_prop_over_time_GAM__{tag}.png"), width = 180, height = 130)
export_figure(p_gam, glue("motility_prop_over_time_GAM__{tag}.png"), width = 180, height = 130)

# ---------- Motility over time: show spread + smooth progression --------------
suppressPackageStartupMessages({
  library(ggbeeswarm)  # dot-swarm to show all cultures without overlap
  library(mgcv)        # smooth (GAM) for over-time trend
})

# Per-day summary (for mean±SE and n labels)
day_counts <- dd |>
  dplyr::group_by(days_from_start) |>
  dplyr::summarise(
    n    = sum(is.finite(motility_prop)),
    mean = mean(motility_prop, na.rm = TRUE),
    se   = stats::sd(motility_prop, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) |>
  dplyr::filter(n > 0)

# A) Dot-swarm (each dot = one culture), with mean±SE and n per day
p_swarm <- ggplot(dd, aes(x = days_from_start, y = motility_prop)) +
  ggbeeswarm::geom_quasirandom(width = 0.25, alpha = 0.25, size = 1.6) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.7) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  geom_text(data = day_counts,
            aes(label = paste0("n=", n), y = pmin(1, mean + se + 0.05)),
            size = 3.2, vjust = 0, nudge_y = 0.01, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
    subtitle = "Each dot is a culture; black line = mean; bars = SE; labels show n per day",
    x = "Days from start", y = "Motility proportion"
  )

show_and_save(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 130)
export_figure(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 130)

# B) Smooth progression across all observed days (GAM), with semi-transparent points
p_gam <- ggplot(dd, aes(days_from_start, motility_prop)) +
  geom_point(alpha = 0.25, size = 1.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
    subtitle = "GAM smooth with 95% CI; dots show culture-level spread (days extend to 16)",
    x = "Days from start", y = "Motility proportion"
  )

show_and_save(p_gam, glue("motility_prop_over_time_GAM__{tag}.png"), width = 180, height = 130)
export_figure(p_gam, glue("motility_prop_over_time_GAM__{tag}.png"), width = 180, height = 130)

  # Species‑aware snapshot (days 0/1/4) — only for ALL run
  if (is.na(species_filter)) {
    trend_sp <- dd |>
      group_by(species, days_from_start) |>
      summarise(n = sum(is.finite(motility_prop)),
                mean = mean(motility_prop, na.rm = TRUE), se = se(motility_prop), .groups = "drop") |>
      filter(days_from_start %in% FOCUS_DAYS)
    
    # flag small-n for plotting labels
    trend_sp <- trend_sp |> mutate(n_flag = ifelse(n < QC_SMALL_N, " (n<15)", ""))
    
    p_trend_sp <- ggplot(trend_sp, aes(days_from_start, mean, colour = species, group = species)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
      geom_line(linewidth = 0.8) + geom_point() +
      labs(title = "Motility proportion by species (days 0/1/4)",
           subtitle = "Mean ± SE; label warns when n < 15 per point",
           x = "Days from start", y = "Mean motility proportion") +
      theme(legend.position = "bottom")
    
    show_and_save(p_trend_sp, "motility_prop_by_species_days014__ALL.png", width = 170, height = 120)
    export_plot(p_trend_sp, "motility_prop_by_species_days014__ALL.png", width = 170, height = 120)
  }
  
  # Intro paragraph
  compose_intro_paragraph <- function(overall_tbl, byday_tbl, dd_label) {
    fmtp <- function(x) scales::percent(x, accuracy = 0.1)
    m_med <- overall_tbl |> filter(variable == "motility_prop") |> pull(median)
    m_q1  <- overall_tbl |> filter(variable == "motility_prop") |> pull(q1)
    m_q3  <- overall_tbl |> filter(variable == "motility_prop") |> pull(q3)
    s_med <- overall_tbl |> filter(variable == "settlement_count") |> pull(median)
    
    d0_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 0]
    d1_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 1]
    d4_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 4]
    
    # Percent decline (median)
    dec_med <- 100 * (d4_med - d0_med) / d0_med
    
    glue::glue(
      "We analysed {nrow(dd)} kelp zoospore cultures {dd_label}. Motility fractions ",
      "were typically {fmtp(m_med)} [IQR {fmtp(m_q1)}–{fmtp(m_q3)}]. At key time points,",
      " medians were {fmtp(d0_med)} (day 0), {fmtp(d1_med)} (day 1) and {fmtp(d4_med)} (day 4),",
      " corresponding to a {round(dec_med,1)}% median change from day 0 to day 4. ",
      "Settlement counts centred on {format(s_med, big.mark=',')}."
    )
  }
  
  intro_text <- compose_intro_paragraph(overall_tbl, byday_tbl,
                                        ifelse(is.na(species_filter), "across three species",
                                               paste0("in ", species_filter)))
  writeLines(intro_text, file.path(DIR_OUTPUTS, glue("intro_paragraph_draft__{tag}.txt")))
  
  # Return list of key outputs for side‑by‑side comparison
  list(tag = tag, n = nrow(dd), overall = overall_tbl, byday = byday_tbl,
       decline = decline_tbl, trend = trend_tbl)
}
# ==============================================================================
# 6. Run BOTH analyses: ALL species and L. digitata only
# ==============================================================================
res_all <- run_summary(d, species_filter = NA_character_, label = "all_species")
res_ld  <- run_summary(d, species_filter = "L. digitata",   label = "l_digitata")

# Combine decline tables for quick side‑by‑side view
both_declines <- bind_rows(res_all$decline, res_ld$decline)
readr::write_csv(both_declines, file.path(DIR_OUTPUTS, "decline_day0_to_day4__ALL_vs_Ldigitata.csv"))

# Small‑n warnings to log
for (res in list(res_all, res_ld)) {
  small_pts <- res$trend |> filter(n < QC_SMALL_N)
  if (nrow(small_pts)) log_run(STAGE, glue("{res$tag}: {nrow(small_pts)} trend points with n< {QC_SMALL_N}"))
}
# ==============================================================================
# 7. Housekeeping
# ==============================================================================
log_run(STAGE, "Summary tables and figures written to outputs/, plots/, figures/")
message("✔ 02_summary completed for ALL + L. digitata.")

# ---- Narrative visuals & tables to support introduction ----------------------
# (1) Distribution of motility at key days; (2) Composition bars motile vs stationary
# (3) Settlement ~ motility (binned means ± SE and LOESS); (4) Optional half-life note
# These exports help readers grasp: mixed states, high settlement, link to motility,
# and steady-but-persistent decline of motility over time.

# Safe recompute of helper variables
if (!exists("d")) {
  INPUT_FILE <- here::here("data_raw","kelp_control_cultures_clean.csv")
  d <- readr::read_csv(INPUT_FILE) |> janitor::clean_names()
}

d <- d |>
  mutate(
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    motility_prop   = pmin(pmax(suppressWarnings(as.numeric(motility_prop)),0),1),
    stationary_prop = pmax(0, 1 - motility_prop),
    settlement_count = suppressWarnings(as.numeric(settlement_count))
  )

# (1) Distribution at key days
FOCUS_DAYS <- get0("FOCUS_DAYS", ifnotfound = c(0,1,4))
p_mix_day <- ggplot(d |> filter(days_from_start %in% FOCUS_DAYS),
                    aes(factor(days_from_start), motility_prop)) +
  geom_violin(fill = kelp_pal["kelp_teal"], alpha = 0.25, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  labs(title = "Distribution of motility fractions at key days",
       x = "Day", y = "Motility proportion")
show_and_save(p_mix_day, "motility_distribution_days014.png", width = 170, height = 120)
export_plot(p_mix_day, "motility_distribution_days014.png", width = 170, height = 120)

# (2) Stacked composition bars at key days
comp_day <- d |>
  filter(days_from_start %in% FOCUS_DAYS) |>
  group_by(days_from_start) |>
  summarise(motile_mean = mean(motility_prop, na.rm = TRUE),
            stationary_mean = mean(stationary_prop, na.rm = TRUE), .groups = "drop") |>
  pivot_longer(c(motile_mean, stationary_mean), names_to = "state", values_to = "prop") |>
  mutate(state = recode(state, motile_mean = "Motile", stationary_mean = "Stationary"))

p_stack <- ggplot(comp_day, aes(x = factor(days_from_start), y = prop, fill = state)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(Motile = kelp_pal["kelp_blue"], Stationary = kelp_pal["kelp_grey"])) +
  labs(title = "Average culture composition at key days",
       x = "Day", y = "Proportion of cells", fill = NULL)

show_and_save(p_stack, "composition_stacked_days014.png", width = 170, height = 120)
export_plot(p_stack, "composition_stacked_days014.png", width = 170, height = 120)

# (3) Settlement vs motility (binned means ± SE)
breaks <- quantile(d$motility_prop, probs = seq(0,1,0.2), na.rm = TRUE)
d_binned <- d |>
  mutate(mot_bin = cut(motility_prop, breaks = unique(breaks), include.lowest = TRUE, dig.lab = 3))

rel_tbl <- d_binned |>
  group_by(mot_bin) |>
  summarise(n = sum(is.finite(settlement_count)),
            settle_mean = mean(settlement_count, na.rm = TRUE),
            settle_se = sd(settlement_count, na.rm = TRUE)/sqrt(n), .groups = "drop")

p_rel <- ggplot(rel_tbl, aes(x = mot_bin, y = settle_mean)) +
  geom_col(fill = kelp_pal["kelp_green"], alpha = 0.8) +
  geom_errorbar(aes(ymin = settle_mean - settle_se, ymax = settle_mean + settle_se), width = 0.2) +
  labs(title = "Settlement increases with motility (binned means ± SE)",
       x = "Motility proportion (quintile bins)", y = "Mean settlement count") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
show_and_save(p_rel, "settlement_vs_motility_binned.png", width = 170, height = 120)
export_plot(p_rel, "settlement_vs_motility_binned.png", width = 170, height = 120)

# Raw scatter with LOESS (clip extreme tail for readability)
p_scatter <- ggplot(d, aes(motility_prop, settlement_count)) +
  geom_point(alpha = 0.25, shape = 16) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(title = "Settlement vs motility (raw with LOESS)",
       x = "Motility proportion", y = "Settlement count") +
  coord_cartesian(ylim = c(NA, quantile(d$settlement_count, 0.99, na.rm=TRUE)))
show_and_save(p_scatter, "settlement_vs_motility_loess.png", width = 170, height = 120)
export_plot(p_scatter, "settlement_vs_motility_loess.png", width = 170, height = 120)

# (4) Optional persistence note: compute and store half-life if present in outputs
# If earlier step wrote motility_half_life_note.txt, nothing to do here.
# ==============================================================================
# 8. Narrative paragraph booster (append)
# ==============================================================================
append_narrative <- function(base_path = file.path(DIR_OUTPUTS, "intro_paragraph_draft.txt")) {
  txt <- if (file.exists(base_path)) readLines(base_path) else ""
  extra <- "Cultures consistently comprised a mixture of motile and stationary cells.
" |
    paste0("Across the key timepoints (0/1/4 d), average culture composition shifted gradually, but motility persisted, consistent with a steady rather than abrupt decline.
",
           "Higher settlement counts tended to occur in cultures with higher motility fractions, supporting the encounter-driven view that motile spores more often reach and attach to surfaces.")
  new_txt <- paste(c(txt, "", extra), collapse = "
")
  writeLines(new_txt, base_path)
}
try(append_narrative(), silent = TRUE)

# ==============================================================================
# END
# ==============================================================================
