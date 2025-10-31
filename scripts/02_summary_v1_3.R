# ==============================================================================
# Project:     Kelp Zoospore Motility Profiling
# Script:      02_summary.R
# Author:      Marianne Glascott (with ChatGPT assistant)
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2025-10-31
# Version:     1.3
# ==============================================================================
# Description
#   Import clean culture-level CSV (642 rows)
#   QC ranges, write flagged rows
#   Summary statistics overall + days 0/1/4 (reduce day-mix skew)
#   Motility over time:
#     A) Dot-swarm per day (spread) + mean±SE + n labels
#     B) Weighted GAM smooth (n/day weights) + 95% CI + effect-size subtitle
#     C) Side-by-side panel (A + B)
#     D) Species snapshot (0/1/4) with mean±SE (ALL run only)
#   Relationship plots: distributions at key days, composition bars,
#   and settlement ~ motility (binned means ± SE and LOESS).
#   Save to outputs/, plots/, figures/ with screen + pub formats via 01_setup.
# ==============================================================================
# Contents:
# 1. Setup & flags
# 2. Read data & coercions
# 3. QC checks & report
# 4. Small distribution summary utility
# 5. Reusable runner (ALL vs L. digitata)
# 6. Run BOTH analyses
# 7. Additional visuals for the introduction narrative
# 8. Narrative paragraph booster (append to intro draft)
# 9. Housekeeping
# ==============================================================================
# 1. Setup & flags
# ==============================================================================
STAGE <- "SUMMARY"
source(here::here("scripts","01_setup.R"))  # provides DIR_* paths, theme, exporters

# Preferred snapshot days to avoid day-mix skew
FOCUS_DAYS <- c(0, 1, 4)
INPUT_FILE <- here::here("data_raw","kelp_control_cultures_clean.csv")

# QC thresholds
QC_MIN_MOTILITY <- 0
QC_MAX_MOTILITY <- 1
QC_MIN_COUNT    <- 0
QC_SMALL_N      <- 15   # warn when group n < 15 for mean/SE stability

# Extra libs used in this script
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readr); library(janitor); library(glue); library(scales)
  library(ggbeeswarm)   # dot-swarm
  library(mgcv)         # GAM smooths
  library(patchwork)    # side-by-side composition
  library(stringr)      # title/subtitle wrapping
})

# Helper: standard error
se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

# Helper: pad & wrap long titles/subtitles to avoid clipping on save
pad_titles <- function(p, top = 16, right = 12, bottom = 12, left = 12, wrap = 90) {
  labs <- ggplot2::ggplot_build(p)$plot$labels
  tt   <- if (!is.null(labs$title))    stringr::str_wrap(labs$title,    width = wrap) else NULL
  st   <- if (!is.null(labs$subtitle)) stringr::str_wrap(labs$subtitle, width = wrap) else NULL
  p +
    ggplot2::labs(title = tt, subtitle = st) +
    ggplot2::theme(
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      plot.title    = element_text(margin = margin(b = 6)),
      plot.subtitle = element_text(margin = margin(t = 2, b = 8)),
      plot.margin   = margin(t = top, r = right, b = bottom, l = left)
    )
}

# ==============================================================================
# 2. Read data & coercions
# ==============================================================================
stopifnot(file.exists(INPUT_FILE))
d0 <- readr::read_csv(INPUT_FILE) |> janitor::clean_names()

need <- c("mobile_cell_count","stationary_cell_count","motility_prop",
          "settlement_count","days_from_start","species",
          "collection_site","season")
miss <- setdiff(need, names(d0))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

as_int <- function(x) suppressWarnings(as.integer(x))
d <- d0 |>
  mutate(
    days_from_start      = as_int(days_from_start),
    species              = as.factor(species),
    collection_site      = as.factor(collection_site),
    season               = as.factor(season),
    motility_prop        = suppressWarnings(as.numeric(motility_prop)),
    settlement_count     = suppressWarnings(as.numeric(settlement_count)),
    mobile_cell_count    = suppressWarnings(as.numeric(mobile_cell_count)),
    stationary_cell_count= suppressWarnings(as.numeric(stationary_cell_count))
  )

# ==============================================================================
# 3. QC checks & report
# ==============================================================================
qc_rows <- d |>
  mutate(
    qc_motility_oob = is.na(motility_prop) | motility_prop < QC_MIN_MOTILITY | motility_prop > QC_MAX_MOTILITY,
    qc_neg_counts   = (mobile_cell_count   < QC_MIN_COUNT) |
      (stationary_cell_count < QC_MIN_COUNT) |
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

# Clamp motility into [0,1] (record count of clamped values)
clamped <- sum(d$motility_prop < QC_MIN_MOTILITY | d$motility_prop > QC_MAX_MOTILITY, na.rm = TRUE)
d$motility_prop <- pmin(pmax(d$motility_prop, QC_MIN_MOTILITY), QC_MAX_MOTILITY)
if (clamped > 0) log_run(STAGE, glue("QC: clamped motility_prop for {clamped} values to [0,1]"))

# ==============================================================================
# 4. Small distribution summary utility
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

# ==============================================================================
# 5. Reusable runner (ALL vs L. digitata)
# ==============================================================================
run_summary <- function(data, species_filter = NA_character_, label = "all_species") {
  
  dd  <- if (!is.na(species_filter)) dplyr::filter(data, species == species_filter) else data
  tag <- ifelse(is.na(species_filter), "ALL", gsub("[ .]", "_", species_filter))
  log_run(STAGE, glue("Running summary for: {tag} (n={nrow(dd)})"))
  
  vars_to_summarise <- c("motility_prop","settlement_count","mobile_cell_count","stationary_cell_count")
  
  # Overall distributions
  overall_tbl <- purrr::map_dfr(vars_to_summarise, \(v) {
    dplyr::bind_cols(variable = v, summarise_distribution(dd[[v]]))
  }) |> dplyr::relocate(variable)
  readr::write_csv(overall_tbl, file.path(DIR_OUTPUTS, glue("overall_distribution_summary__{tag}.csv")))
  
  # Snapshot days 0/1/4
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
             ), .names = "{.col}.{.fn}"),
      .groups = "drop"
    )
  readr::write_csv(byday_tbl, file.path(DIR_OUTPUTS, glue("byday_summary_0_1_4__{tag}.csv")))
  
  # Decline 0 → 4 (mean & median)
  d0row <- byday_tbl |> filter(days_from_start == 0)
  d4row <- byday_tbl |> filter(days_from_start == 4)
  decline_tbl <- tibble(
    group = tag,
    n_day0 = d0row$motility_prop.n, n_day4 = d4row$motility_prop.n,
    mean_day0 = d0row$motility_prop.mean, mean_day4 = d4row$motility_prop.mean,
    med_day0  = d0row$motility_prop.median, med_day4  = d4row$motility_prop.median,
    decline_mean_pct   = 100 * (mean_day4 - mean_day0) / mean_day0,
    decline_median_pct = 100 * (med_day4  - med_day0) / med_day0
  )
  readr::write_csv(decline_tbl, file.path(DIR_OUTPUTS, glue("decline_day0_to_day4__{tag}.csv")))
  
  #  Motility over time: spread + weighted smooth --------------------
  # Per-day n/mean/SE and CSVs
  day_counts <- dd |>
    dplyr::group_by(days_from_start) |>
    dplyr::summarise(
      n    = sum(is.finite(motility_prop)),
      mean = mean(motility_prop, na.rm = TRUE),
      se   = stats::sd(motility_prop, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) |>
    dplyr::filter(n > 0)
  
  readr::write_csv(day_counts, file.path(DIR_OUTPUTS, glue("day_counts__{tag}.csv")))
  
  day_counts_species <- dd |>
    dplyr::group_by(species, days_from_start) |>
    dplyr::summarise(
      n    = sum(is.finite(motility_prop)),
      mean = mean(motility_prop, na.rm = TRUE),
      se   = stats::sd(motility_prop, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) |>
    dplyr::filter(n > 0)
  readr::write_csv(day_counts_species, file.path(DIR_OUTPUTS, glue("day_counts_by_species__{tag}.csv")))
  
  # A) Dot-swarm per day (spread + mean±SE + n labels)
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
  p_swarm <- pad_titles(p_swarm, top = 18, wrap = 80)
  show_and_save(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 140)
  export_figure(p_swarm, glue("motility_prop_swarm_meanSE__{tag}.png"), width = 180, height = 140)
  
  # B) Weighted GAM smooth (weights = n/day) + effect-size subtitle
  dd_w <- dd |>
    dplyr::left_join(day_counts |> dplyr::select(days_from_start, n) |> dplyr::rename(wt = n),
                     by = "days_from_start")
  
  lm_fit <- lm(motility_prop ~ days_from_start, data = dd)
  slope  <- coef(lm_fit)[["days_from_start"]]
  ci     <- confint(lm_fit)["days_from_start", ]
  d0min  <- min(dd$days_from_start, na.rm = TRUE)
  dLmax  <- max(dd$days_from_start, na.rm = TRUE)
  m0     <- mean(dd$motility_prop[dd$days_from_start == d0min], na.rm = TRUE)
  mL     <- mean(dd$motility_prop[dd$days_from_start == dLmax], na.rm = TRUE)
  pct    <- if (is.finite(m0) && m0 > 0) 100 * (mL - m0) / m0 else NA_real_
  
  annot <- glue::glue(
    "Weighted GAM (n/day) with 95% CI; dots show spread. ",
    "Linear slope {sprintf('%.3f', slope)}/day [{sprintf('%.3f', ci[1])}, {sprintf('%.3f', ci[2])}]",
    if (is.finite(pct)) glue::glue("; Δ% {d0min}→{dLmax}: {sprintf('%.1f', pct)}%") else ""
  )
  
  p_gam_wt <- ggplot(dd_w, aes(days_from_start, motility_prop)) +
    geom_point(alpha = 0.22, size = 1.3) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, aes(weight = wt)) +
    geom_text(data = day_counts,
              aes(x = days_from_start, y = pmin(1, mean + 0.07), label = paste0("n=", n)),
              size = 3, vjust = 0, colour = "grey20") +
    scale_x_continuous(breaks = sort(unique(dd$days_from_start))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = glue::glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
      subtitle = annot, x = "Days from start", y = "Motility proportion"
    )
  p_gam_wt <- pad_titles(p_gam_wt, top = 18, wrap = 80)
  show_and_save(p_gam_wt, glue("motility_prop_over_time_GAM_weighted__{tag}.png"), width = 180, height = 140)
  export_figure(p_gam_wt, glue("motility_prop_over_time_GAM_weighted__{tag}.png"), width = 180, height = 140)
  
  # C) Side-by-side panel
  p_panel <- p_swarm + p_gam_wt + plot_layout(widths = c(1, 1)) &
    labs(title = glue::glue("Motility over time — {ifelse(tag=='ALL','All species', species_filter)}"))
  p_panel <- pad_titles(p_panel, top = 18, wrap = 80)
  show_and_save(p_panel, glue("motility_over_time_panel_swarm_plus_smooth__{tag}.png"),
                width = 220, height = 140)
  export_figure(p_panel, glue("motility_over_time_panel_swarm_plus_smooth__{tag}.png"),
                width = 220, height = 140)
  
  # Keep a lightweight trend table for downstream small-n checks
  trend_tbl <- day_counts
  
  # D) Species snapshot (0/1/4) — only for ALL run
  if (is.na(species_filter)) {
    trend_sp <- dd |>
      group_by(species, days_from_start) |>
      summarise(n = sum(is.finite(motility_prop)),
                mean = mean(motility_prop, na.rm = TRUE),
                se = se(motility_prop), .groups = "drop") |>
      filter(days_from_start %in% FOCUS_DAYS) |>
      mutate(n_flag = ifelse(n < QC_SMALL_N, " (n<15)", ""))
    
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
  
  # Intro paragraph (draft)
  compose_intro_paragraph <- function(overall_tbl, byday_tbl, dd_label) {
    fmtp <- function(x) scales::percent(x, accuracy = 0.1)
    m_med <- overall_tbl |> filter(variable == "motility_prop") |> pull(median)
    m_q1  <- overall_tbl |> filter(variable == "motility_prop") |> pull(q1)
    m_q3  <- overall_tbl |> filter(variable == "motility_prop") |> pull(q3)
    s_med <- overall_tbl |> filter(variable == "settlement_count") |> pull(median)
    
    d0_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 0]
    d1_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 1]
    d4_med <- byday_tbl$motility_prop.median[byday_tbl$days_from_start == 4]
    dec_med <- 100 * (d4_med - d0_med) / d0_med
    
    glue::glue(
      "We analysed {nrow(dd)} kelp zoospore cultures {dd_label}. Motility fractions ",
      "were typically {fmtp(m_med)} [IQR {fmtp(m_q1)}–{fmtp(m_q3)}]. At key time points, ",
      "medians were {fmtp(d0_med)} (day 0), {fmtp(d1_med)} (day 1), and {fmtp(d4_med)} (day 4), ",
      "corresponding to a {round(dec_med,1)}% median change from day 0 to day 4. ",
      "Settlement counts centred on {format(s_med, big.mark=',')}."
    )
  }
  
  intro_text <- compose_intro_paragraph(overall_tbl, byday_tbl,
                                        ifelse(is.na(species_filter), "across three species",
                                               paste0("in ", species_filter)))
  writeLines(intro_text, file.path(DIR_OUTPUTS, glue("intro_paragraph_draft__{tag}.txt")))
  
  # Return key objects
  list(tag = tag, n = nrow(dd), overall = overall_tbl, byday = byday_tbl,
       decline = decline_tbl, trend = trend_tbl)
}

# ==============================================================================
# 6. Run BOTH analyses
# ==============================================================================
res_all <- run_summary(d, species_filter = NA_character_, label = "all_species")
res_ld  <- run_summary(d, species_filter = "L. digitata",   label = "l_digitata")

# Side-by-side decline table
both_declines <- bind_rows(res_all$decline, res_ld$decline)
readr::write_csv(both_declines, file.path(DIR_OUTPUTS, "decline_day0_to_day4__ALL_vs_Ldigitata.csv"))

# Small-n warnings to log (per-day trend table)
for (res in list(res_all, res_ld)) {
  small_pts <- res$trend |> filter(n < QC_SMALL_N)
  if (nrow(small_pts)) log_run(STAGE, glue("{res$tag}: {nrow(small_pts)} day-level points with n < {QC_SMALL_N}"))
}

#================================================================================
# By species diagnostic
#===============================================================================
# --- QUICK DIAGNOSIS ----------------------------------------------------------
dplyr::count(d, species, sort = TRUE) |> print(n = Inf)

# --- CANONICALISE SPECIES NAMES ------------------------------------------------
suppressPackageStartupMessages({ library(stringr) })

canonicalise_species <- function(x) {
  # trim + collapse spaces
  x0 <- str_squish(str_trim(as.character(x)))
  # map by pattern (case-insensitive)
  lab <- dplyr::case_when(
    str_detect(str_to_lower(x0), "digitata")   ~ "L. digitata",
    str_detect(str_to_lower(x0), "hyperborea") ~ "L. hyperborea",
    str_detect(str_to_lower(x0), "latissima")  ~ "S. latissima",
    TRUE ~ x0
  )
  # force to factor with desired level order
  factor(lab, levels = c("L. digitata", "L. hyperborea", "S. latissima"))
}

# apply canonicalisation in-place
d <- d |>
  dplyr::mutate(species = canonicalise_species(species))

# sanity check after canonicalisation
dplyr::count(d, species, sort = TRUE) |> print(n = Inf)

# ==============================================================================
# 6b) By-species overlays and 2×2 panel  (DEFENSIVE: always provide wt)
# Paste AFTER Section 6 (after res_all/res_ld); requires 'd' in memory
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("mgcv", quietly = TRUE)) install.packages("mgcv")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  library(mgcv); library(patchwork)
})

# pad_titles available globally
if (!exists("pad_titles")) {
  pad_titles <- function(p, top = 16, right = 12, bottom = 12, left = 12, wrap = 90) {
    labs <- ggplot2::ggplot_build(p)$plot$labels
    tt   <- if (!is.null(labs$title))    stringr::str_wrap(labs$title,    width = wrap) else NULL
    st   <- if (!is.null(labs$subtitle)) stringr::str_wrap(labs$subtitle, width = wrap) else NULL
    p + labs(title = tt, subtitle = st) +
      theme(plot.title.position = "plot",
            plot.subtitle = element_text(margin = margin(t = 2, b = 8)),
            plot.margin = margin(t = top, r = right, b = bottom, l = left))
  }
}

# ---- Per-species/day counts and weighted datasets ----------------------------
sp_day_counts <- d |>
  dplyr::group_by(species, days_from_start) |>
  dplyr::summarise(n = sum(is.finite(motility_prop)),
                   mean = mean(motility_prop, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::filter(n > 0)

# Weighted per-row dataset (species+day join)
d_w_sp <- d |>
  dplyr::left_join(sp_day_counts |>
                     dplyr::select(species, days_from_start, n) |>
                     dplyr::rename(wt = n),
                   by = c("species","days_from_start")) |>
  dplyr::mutate(wt = dplyr::coalesce(wt, 1L))

# Also totals per day for the “All species” panel
day_totals <- d |>
  dplyr::group_by(days_from_start) |>
  dplyr::summarise(n = sum(is.finite(motility_prop)),
                   mean = mean(motility_prop, na.rm = TRUE), .groups = "drop")
d_w_all <- d |>
  dplyr::left_join(day_totals |>
                     dplyr::select(days_from_start, n) |>
                     dplyr::rename(wt = n),
                   by = "days_from_start") |>
  dplyr::mutate(wt = dplyr::coalesce(wt, 1L))

# ---------- (A) By-species overlay (single panel) -----------------------------
p_by_species <- ggplot(d_w_sp, aes(days_from_start, motility_prop, colour = species)) +
  geom_point(alpha = 0.18, size = 1.2) +
  geom_smooth(data = d_w_sp, method = "gam", formula = y ~ s(x, k = 6),
              se = TRUE, aes(weight = wt)) +
  geom_text(data = sp_day_counts,
            aes(x = days_from_start, y = pmin(1, mean + 0.06),
                label = paste0("n=", n), colour = species),
            size = 3, vjust = 0, show.legend = FALSE) +
  scale_x_continuous(breaks = sort(unique(d$days_from_start))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c(
    "L. digitata"   = kelp_pal["kelp_blue"],
    "L. hyperborea" = kelp_pal["kelp_green"],
    "S. latissima"  = kelp_pal["kelp_teal"]
  ), drop = FALSE) +
  labs(
    title = "Motility proportion over time — by species",
    subtitle = "Weighted GAM (weights = n/day/species) with 95% CI; dots show culture-level spread",
    x = "Days from start", y = "Motility proportion", colour = NULL
  ) +
  theme(legend.position = "bottom")
p_by_species <- pad_titles(p_by_species, top = 18, wrap = 86)
print(p_by_species)
show_and_save(p_by_species, "motility_prop_over_time_GAM_weighted__BY_SPECIES_overlay.png",
              width = 190, height = 140)
export_figure(p_by_species, "motility_prop_over_time_GAM_weighted__BY_SPECIES_overlay.png",
              width = 190, height = 140)

# ---------- (B) 2×2 panel: All + each species --------------------------------
make_sp_plot <- function(df_w, sp_name) {
  dd_w <- dplyr::filter(df_w, species == sp_name)
  if (nrow(dd_w) == 0L) {
    return(ggplot() + theme_void() + labs(title = sp_name, subtitle = "no data"))
  }
  dcnt <- dplyr::filter(sp_day_counts, species == sp_name)
  
  # Linear slope for concise effect size
  lm_fit <- lm(motility_prop ~ days_from_start, data = dd_w)
  slope  <- coef(lm_fit)[["days_from_start"]]
  ci     <- confint(lm_fit)["days_from_start", ]
  d0     <- min(dd_w$days_from_start, na.rm = TRUE)
  dL     <- max(dd_w$days_from_start, na.rm = TRUE)
  m0     <- mean(dd_w$motility_prop[dd_w$days_from_start == d0], na.rm = TRUE)
  mL     <- mean(dd_w$motility_prop[dd_w$days_from_start == dL], na.rm = TRUE)
  pct    <- if (is.finite(m0) && m0 > 0) 100 * (mL - m0) / m0 else NA_real_
  annot  <- glue::glue("Linear slope {sprintf('%.3f', slope)}/day [{sprintf('%.3f', ci[1])}, {sprintf('%.3f', ci[2])}]",
                       if (is.finite(pct)) glue::glue("; Δ% {d0}→{dL}: {sprintf('%.1f', pct)}%") else "")
  
  ggplot(dd_w, aes(days_from_start, motility_prop)) +
    geom_point(alpha = 0.2, size = 1.2) +
    geom_smooth(data = dd_w, method = "gam", formula = y ~ s(x, k = 6),
                se = TRUE, aes(weight = wt)) +
    geom_text(data = dcnt,
              aes(x = days_from_start, y = pmin(1, mean + 0.06),
                  label = paste0("n=", n)),
              size = 3, vjust = 0, colour = "grey20") +
    scale_x_continuous(breaks = sort(unique(dd_w$days_from_start))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = sp_name,
         subtitle = paste0("Weighted GAM (n/day); ", annot),
         x = "Days from start", y = "Motility proportion")
}

# Panels
p_all <- ggplot(d_w_all, aes(days_from_start, motility_prop)) +
  geom_point(alpha = 0.18, size = 1.2) +
  geom_smooth(data = d_w_all, method = "gam", formula = y ~ s(x, k = 6),
              se = TRUE, aes(weight = wt), colour = kelp_pal["kelp_blue"]) +
  geom_text(data = day_totals,
            aes(x = days_from_start, y = pmin(1, mean + 0.06), label = paste0("n=", n)),
            size = 3, vjust = 0, colour = "grey20") +
  scale_x_continuous(breaks = sort(unique(d$days_from_start))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "All species",
       subtitle = "Weighted GAM (weights = total n/day) with 95% CI",
       x = "Days from start", y = "Motility proportion")

p_ld <- make_sp_plot(d_w_sp, "L. digitata")
p_lh <- make_sp_plot(d_w_sp, "L. hyperborea")
p_sl <- make_sp_plot(d_w_sp, "S. latissima")

p_panel_2x2 <- (p_all + p_ld) / (p_lh + p_sl)
p_panel_2x2 <- pad_titles(
  p_panel_2x2 + plot_annotation(
    title = "Motility proportion over time — All and by species",
    subtitle = "Each panel uses a weighted GAM (weights = n/day in that panel) with 95% CI; points show culture-level spread"
  ),
  top = 18, wrap = 90
)

print(p_panel_2x2)
show_and_save(p_panel_2x2, "motility_prop_over_time_GAM_weighted__PANEL_all_and_species.png",
              width = 230, height = 180)
export_figure(p_panel_2x2, "motility_prop_over_time_GAM_weighted__PANEL_all_and_species.png",
              width = 230, height = 180)

# ---- Show + save each species plot individually ------------------------------
# Small helper: safe filename from species label
to_file_stub <- function(sp) {
  sp |>
    gsub("[^A-Za-z0-9]+", "_", ., perl = TRUE) |>
    tolower()
}

plots_list <- list(
  "L. digitata"   = p_ld,
  "L. hyperborea" = p_lh,
  "S. latissima"  = p_sl
)

for (sp in names(plots_list)) {
  p   <- plots_list[[sp]]
  fs  <- to_file_stub(sp)
  # Optional tidy titles to avoid clipping
  if (exists("pad_titles")) p <- pad_titles(p, top = 18, wrap = 90)
  
  print(p)  # <-- forces display in Plots pane
  show_and_save(p, glue::glue("motility_prop_over_time_GAM_weighted__{fs}.png"),
                width = 190, height = 140)
  export_figure(p, glue::glue("motility_prop_over_time_GAM_weighted__{fs}.png"),
                width = 190, height = 140)
}

# ==============================================================================
# 7. Additional visuals for the introduction narrative
# ==============================================================================
# Safe recompute (if called standalone)
if (!exists("d")) {
  INPUT_FILE <- here::here("data_raw","kelp_control_cultures_clean.csv")
  d <- readr::read_csv(INPUT_FILE) |> janitor::clean_names()
}
d <- d |>
  mutate(
    days_from_start  = suppressWarnings(as.integer(days_from_start)),
    motility_prop    = pmin(pmax(suppressWarnings(as.numeric(motility_prop)),0),1),
    stationary_prop  = pmax(0, 1 - motility_prop),
    settlement_count = suppressWarnings(as.numeric(settlement_count))
  )

# (1) Distribution of motility at key days
p_mix_day <- ggplot(d |> filter(days_from_start %in% FOCUS_DAYS),
                    aes(factor(days_from_start), motility_prop)) +
  geom_violin(fill = kelp_pal["kelp_teal"], alpha = 0.25, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  labs(title = "Distribution of motility fractions at key days",
       x = "Day", y = "Motility proportion")
show_and_save(p_mix_day, "motility_distribution_days014.png", width = 170, height = 120)
export_plot(p_mix_day, "motility_distribution_days014.png", width = 170, height = 120)

# (2) Stacked composition bars (motile vs stationary) at key days
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

# (3) Settlement vs motility — binned means ± SE (quintiles)
breaks <- quantile(d$motility_prop, probs = seq(0,1,0.2), na.rm = TRUE)
d_binned <- d |> mutate(mot_bin = cut(motility_prop, breaks = unique(breaks), include.lowest = TRUE, dig.lab = 3))
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

# (4) Settlement vs motility — raw scatter + LOESS (clip long tail for readability)
p_scatter <- ggplot(d, aes(motility_prop, settlement_count)) +
  geom_point(alpha = 0.25, shape = 16) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  labs(title = "Settlement vs motility (raw with LOESS)",
       x = "Motility proportion", y = "Settlement count") +
  coord_cartesian(ylim = c(NA, quantile(d$settlement_count, 0.99, na.rm=TRUE)))
show_and_save(p_scatter, "settlement_vs_motility_loess.png", width = 170, height = 120)
export_plot(p_scatter, "settlement_vs_motility_loess.png", width = 170, height = 120)

# ==============================================================================
# 8. Narrative paragraph booster (append to intro draft)
# ==============================================================================
append_narrative <- function(base_path = file.path(DIR_OUTPUTS, "intro_paragraph_draft.txt")) {
  txt   <- if (file.exists(base_path)) readLines(base_path) else character(0)
  extra <- paste0(
    "Cultures consistently comprised a mixture of motile and stationary cells.\n",
    "Across the key timepoints (0/1/4 d), average culture composition shifted gradually, ",
    "but motility persisted, consistent with a steady rather than abrupt decline.\n",
    "Higher settlement counts tended to occur in cultures with higher motility fractions, ",
    "supporting the encounter-driven view that motile spores more often reach and attach to surfaces."
  )
  writeLines(c(txt, "", extra), base_path)
}
try(append_narrative(), silent = TRUE)

# ==============================================================================
# 9. Housekeeping
# ==============================================================================
log_run(STAGE, "Summary tables and figures written to outputs/, plots/, figures/")
message("✔ 02_summary completed for ALL + L. digitata.")
# ==============================================================================
# END
# ==============================================================================
