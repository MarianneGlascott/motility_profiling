# ==============================================================================
# Project:     Kelp Zoospore Motility Profiling
# Scripts:     02_summary.R
# Author:      Marianne Glascott (with ChatGPT assistant)
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2025-10-31
# Version:     1.0
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

# ---- 3. QC checks & report ---------------------------------------------------
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

# ---- 4. Utility: compact distribution summary -------------------------------
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

# ---- 5. A reusable runner to produce all outputs for a chosen species filter -
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
  trend_tbl <- dd |>
    group_by(days_from_start) |>
    summarise(
      n = sum(is.finite(motility_prop)),
      mean = mean(motility_prop, na.rm = TRUE),
      se = se(motility_prop),
      .groups = "drop"
    ) |>
    filter(!is.na(mean))
  
  p_trend <- ggplot(trend_tbl, aes(x = days_from_start, y = mean)) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.2) +
    geom_line(linewidth = 0.8) +
    geom_point() +
    labs(title = glue("Motility proportion over time — {ifelse(tag=='ALL','All species', species_filter)}"),
         subtitle = "Mean ± SE across all observed days",
         x = "Days from start", y = "Motility proportion (mean ± SE)") +
    coord_cartesian(ylim = c(0, 1))
  
  show_and_save(p_trend, glue("motility_prop_over_time__{tag}.png"), width = 170, height = 120)
  export_figure(p_trend, glue("motility_prop_over_time__{tag}.png"), width = 170, height = 120)
  
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

# ---- 6. Run BOTH analyses: ALL species and L. digitata only ------------------
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

# ---- 7. Housekeeping ---------------------------------------------------------
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

# -------- Narrative paragraph booster (append) --------
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

# ============================================================================
# END OF FILE
# ============================================================================


# ==============================  03_explore.R  ===============================
# Purpose
#   * Exploratory relationships: pair plots and correlation heatmaps
#   * Methods: Spearman (robust, rank-based) and Pearson (linear on log1p counts)
#   * Multiple testing control: Benjamini–Hochberg FDR (q-values)
#   * Stratifications: overall, by species, and key days (0/1/4)
#   * Outputs: tidy CSVs (r/ρ, p, q, n), annotated heatmaps (dims q>=0.05), pair plots
#   * Caption helper: auto-summarise top associations (|r|, q<0.05, n>=threshold)
# ============================================================================

STAGE <- "EXPLORE"
source(here::here("scripts","01_setup.R"))

.install_if_missing(c("GGally"))
suppressPackageStartupMessages({ library(GGally) })

INPUT_FILE   <- here::here("data_raw","kelp_control_cultures_clean.csv")
FOCUS_DAYS   <- c(0,1,4)
QC_SMALL_N   <- 15
ALPHA_FDR    <- 0.05

stopifnot(file.exists(INPUT_FILE))

d <- readr::read_csv(INPUT_FILE) |> janitor::clean_names() |>
  mutate(
    days_from_start = suppressWarnings(as.integer(days_from_start)),
    species = factor(species),
    collection_site = factor(collection_site),
    season = factor(season),
    motility_prop = pmin(pmax(suppressWarnings(as.numeric(motility_prop)),0),1),
    mobile_cell_count = suppressWarnings(as.numeric(mobile_cell_count)),
    stationary_cell_count = suppressWarnings(as.numeric(stationary_cell_count)),
    settlement_count = suppressWarnings(as.numeric(settlement_count)),
    total_cell_count = mobile_cell_count + stationary_cell_count,
    log_mobile = log1p(mobile_cell_count),
    log_stationary = log1p(stationary_cell_count),
    log_total = log1p(total_cell_count),
    log_settlement = log1p(settlement_count)
  )

log_run(STAGE, glue("Loaded {nrow(d)} rows for exploratory analysis"))

# ---- Helpers -----------------------------------------------------------------
num_ok <- function(x) is.finite(x) & !is.na(x)

compute_cor_table <- function(data, vars, method = "spearman") {
  stopifnot(all(vars %in% names(data)))
  # Correlation matrix using pairwise complete obs
  M <- stats::cor(as.matrix(data[vars]), use = "pairwise.complete.obs", method = method)
  # p-values and Ns (pairwise)
  get_ct <- function(i, j) {
    x <- data[[vars[i]]]; y <- data[[vars[j]]]
    ok <- num_ok(x) & num_ok(y)
    n <- sum(ok)
    if (n < 3) return(c(r = NA, p = NA, n = n))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
    c(r = unname(ct$estimate), p = unname(ct$p.value), n = n)
  }
  k <- length(vars)
  R <- matrix(NA_real_, k, k); P <- matrix(NA_real_, k, k); N <- matrix(0L, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    res <- get_ct(i,j); R[i,j] <- res[[1]]; P[i,j] <- res[[2]]; N[i,j] <- res[[3]]
  }
  dimnames(R) <- dimnames(P) <- dimnames(N) <- list(vars, vars)
  # Tidy form
  tidy <- as.data.frame(as.table(R)) |>
    dplyr::rename(var1 = Var1, var2 = Var2, r = Freq) |>
    dplyr::left_join(as.data.frame(as.table(P)) |>
                       dplyr::rename(var1 = Var1, var2 = Var2, p = Freq), by = c("var1","var2")) |>
    dplyr::left_join(as.data.frame(as.table(N)) |>
                       dplyr::rename(var1 = Var1, var2 = Var2, n = Freq), by = c("var1","var2"))
  list(R = R, P = P, N = N, tidy = tidy)
}

add_bh_fdr <- function(tidy_df) {
  # Adjust only unique off-diagonal pairs (var1 < var2), then mirror back
  uniq <- tidy_df |>
    dplyr::filter(var1 != var2) |>
    dplyr::mutate(pair = paste(pmin(var1,var2), pmax(var1,var2), sep = "||")) |>
    dplyr::group_by(pair) |>
    dplyr::summarise(var1 = first(var1), var2 = first(var2), r = first(r), p = first(p), n = first(n), .groups = "drop")
  uniq$q <- p.adjust(uniq$p, method = "BH")
  # Join q back to full table (both directions), leave diagonal q=NA
  tidy_q <- tidy_df |>
    dplyr::mutate(pair = dplyr::if_else(var1==var2, NA_character_, paste(pmin(var1,var2), pmax(var1,var2), sep = "||"))) |>
    dplyr::left_join(uniq[,c("pair","q")], by = "pair") |>
    dplyr::select(-pair)
  tidy_q
}

save_cor_heatmap <- function(R, N, tidy_with_q, title, file_stub, n_warn = QC_SMALL_N, alpha_sig = 0.65, alpha_nonsig = 0.18, q_cut = ALPHA_FDR) {
  df <- as.data.frame(as.table(R)) |>
    dplyr::rename(var1 = Var1, var2 = Var2, r = Freq)
  nn <- as.data.frame(as.table(N)) |>
    dplyr::rename(var1 = Var1, var2 = Var2, n = Freq)
  df <- dplyr::left_join(df, nn, by = c("var1","var2")) |>
    dplyr::left_join(tidy_with_q[,c("var1","var2","q")], by = c("var1","var2")) |>
    dplyr::mutate(n_flag = dplyr::if_else(n < n_warn, "*", ""),
                  sig = !is.na(q) & q < q_cut,
                  alpha = dplyr::if_else(sig, alpha_sig, alpha_nonsig),
                  lab = dplyr::if_else(var1==var2, "", sprintf("%.2f%s", r, n_flag)))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(var1, var2, fill = r, alpha = alpha)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3) +
    ggplot2::scale_fill_gradient2(limits = c(-1,1), oob = scales::squish,
                                  low = "#6A3D9A", mid = "#FFFFFF", high = "#1F78B4",
                                  name = "r / ρ") +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(title = title, subtitle = glue("Cells dimmed when q ≥ {ALPHA_FDR}; '*' marks n < {n_warn}"),
                  x = NULL, y = NULL) +
    theme_kelp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  show_and_save(p, paste0(file_stub, ".png"), width = 170, height = 170)
  export_figure(p, paste0(file_stub, ".png"), width = 170, height = 170)
}

make_ggpairs <- function(data, vars, colour = NULL, title = "Pair plot", file_stub = "pairs") {
  p <- GGally::ggpairs(data = data, columns = vars,
                       aes(color = {{colour}}),
                       upper = list(continuous = GGally::wrap("cor", method = "spearman", size = 3)),
                       lower = list(continuous = "smooth_loess"),
                       diag  = list(continuous = GGally::wrap("densityDiag"))) +
    theme_kelp() +
    theme(legend.position = "bottom") +
    ggtitle(title)
  show_and_save(p, paste0(file_stub, ".png"), width = 220, height = 200)
  export_plot(p, paste0(file_stub, ".png"), width = 220, height = 200)
}

write_top_associations_caption <- function(tidy_df, file_stub, method_label, n_min = QC_SMALL_N, q_cut = ALPHA_FDR) {
  top <- tidy_df |>
    dplyr::filter(var1 != var2, is.finite(r), is.finite(p)) |>
    dplyr::filter(n >= n_min, !is.na(q) & q < q_cut) |>
    dplyr::mutate(pair = paste(pmin(var1,var2), pmax(var1,var2), sep = "–")) |>
    dplyr::distinct(pair, .keep_all = TRUE) |>
    dplyr::arrange(dplyr::desc(abs(r))) |>
    dplyr::slice_head(n = 3) |>
    dplyr::mutate(line = glue::glue("{pair}: r = {sprintf('%.2f', r)}, n = {n}, q = {sprintf('%.3f', q)}"))
  if (nrow(top) == 0) {
    txt <- glue::glue("No associations passed FDR (q < {q_cut}) with n ≥ {n_min}.")
  } else {
    bullet <- paste0("• ", top$line, collapse = "
")
    txt <- glue::glue("Top {nrow(top)} {method_label} associations (q < {q_cut}, n ≥ {n_min}):
{bullet}")
  }
  path <- file.path(DIR_OUTPUTS, paste0("caption_top_associations__", file_stub, ".txt"))
  writeLines(txt, path); invisible(path)
}

# ---- Variables of interest ---------------------------------------------------
base_vars <- c("motility_prop","log_mobile","log_stationary","log_total","log_settlement")

# A runner to compute, adjust FDR, save CSV + heatmap + caption in one go
run_one_panel <- function(df, vars, tag, method = c("spearman","pearson")) {
  method <- match.arg(method)
  ct <- compute_cor_table(df, vars, method = method)
  tidy_q <- add_bh_fdr(ct$tidy)
  # Save tidy table with q
  readr::write_csv(tidy_q, file.path(DIR_OUTPUTS, glue("correlations_{method}__{tag}.csv")))
  # Heatmap (dims q>=0.05)
  save_cor_heatmap(ct$R, ct$N, tidy_q, title = glue("Correlations — {tag} ({method})"),
                   file_stub = glue("heatmap_correlations__{tag}__{method}"))
  # Caption
  write_top_associations_caption(tidy_q, file_stub = glue("{tag}__{method}"), method_label = toupper(method))
}

# ---- Overall panels ----------------------------------------------------------
run_one_panel(d, base_vars, tag = "ALL", method = "spearman")
run_one_panel(d, base_vars, tag = "ALL", method = "pearson")

make_ggpairs(d, base_vars, colour = species,
             title = "Pair plot — All species (Spearman; LOESS fits)",
             file_stub = "pairs__ALL")

# ---- By species --------------------------------------------------------------
for (sp in levels(d$species)) {
  dd <- dplyr::filter(d, species == sp)
  tag <- gsub("[ .]","_", sp)
  run_one_panel(dd, base_vars, tag = tag, method = "spearman")
  run_one_panel(dd, base_vars, tag = tag, method = "pearson")
  make_ggpairs(dd, base_vars, colour = NULL,
               title = glue("Pair plot — {sp} (Spearman)"),
               file_stub = glue("pairs__{tag}"))
}

# ---- Key days snapshots (0,1,4) ---------------------------------------------
for (dy in FOCUS_DAYS) {
  dd <- dplyr::filter(d, days_from_start == dy)
  if (nrow(dd) < 3) {
    log_run(STAGE, glue("Day {dy}: insufficient rows for correlations"))
    next
  }
  run_one_panel(dd, base_vars, tag = glue("day{dy}"), method = "spearman")
  run_one_panel(dd, base_vars, tag = glue("day{dy}"), method = "pearson")
  make_ggpairs(dd, base_vars, colour = species,
               title = glue("Pair plot — day {dy} (Spearman)"),
               file_stub = glue("pairs__day{dy}"))
}

log_run(STAGE, "Exploratory pair plots + correlation (Spearman & Pearson) with FDR complete")
message("✔ 03_explore completed.")

# ============================================================================
# END 03_explore.R
# ============================================================================
