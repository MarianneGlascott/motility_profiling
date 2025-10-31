# ==============================================================================
# Project:     Kelp Zoospore Motility Profiling
# Scripts:     03_explore.R
# Author:      Marianne Glascott (with ChatGPT assistant)
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2025-10-31
# Version:     1.0
# ==============================================================================
# Description
#   Exploratory relationships: pair plots and correlation heatmaps
#   Methods: Spearman (robust, rank-based) and Pearson (linear on log1p counts)
#   Multiple testing control: Benjamini–Hochberg FDR (q-values)
#   Stratifications: overall, by species, and key days (0/1/4)
#   Outputs: tidy CSVs (r/ρ, p, q, n), annotated heatmaps (dims q>=0.05), pair plots
#   Caption helper: auto-summarise top associations (|r|, q<0.05, n>=threshold)
# ==============================================================================
# Notes
# - Paths are project‑relative using here::here().
# - Figures saved to plots/ (screen/EDA) and figures/ (pub‑quality) in PNG, TIFF, PDF.
# - Git: points origin to https://github.com/MarianneGlascott/motility_profiling.git
# - Designed to work from RStudio Project: C:/Motility_Profiling/Motility_Profiling.Rproj
# ==============================================================================
# Contents:
# 1. setup
# 2. Helpers
# 3. Variables of interest 
# 4. Overall panels 
# 5. By species 
# 6. Key days snapshots (0,1,4)
# ==============================================================================
# 1. setup
#===============================================================================
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
#===============================================================================
# 2. Helpers
#===============================================================================
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
#===============================================================================
# 3. Variables of interest 
#===============================================================================
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
#===============================================================================
# 4. Overall panels
#===============================================================================
run_one_panel(d, base_vars, tag = "ALL", method = "spearman")
run_one_panel(d, base_vars, tag = "ALL", method = "pearson")

make_ggpairs(d, base_vars, colour = species,
             title = "Pair plot — All species (Spearman; LOESS fits)",
             file_stub = "pairs__ALL")
#===============================================================================
# 5. By species 
#===============================================================================
for (sp in levels(d$species)) {
  dd <- dplyr::filter(d, species == sp)
  tag <- gsub("[ .]","_", sp)
  run_one_panel(dd, base_vars, tag = tag, method = "spearman")
  run_one_panel(dd, base_vars, tag = tag, method = "pearson")
  make_ggpairs(dd, base_vars, colour = NULL,
               title = glue("Pair plot — {sp} (Spearman)"),
               file_stub = glue("pairs__{tag}"))
}
#===============================================================================
# 6. Key days snapshots (0,1,4)
#===============================================================================
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

#===============================================================================
#===============================================================================
#===============================================================================
# Additional plots to build on what we have already
#===============================================================================
# Additional packages & data
#===============================================================================
# PACKAGES

need <- c("mgcv","glmmTMB","dplyr","readr","ppcor")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
lapply(need, library, character.only = TRUE)

## ---- Load controls file (your path) ----
# Use forward slashes on Windows:
dat_path <- "C:/Motility_Profiling/data_raw/kelp_control_cultures_clean.csv"
d_spore  <- readr::read_csv(dat_path, guess_max = 100000)

## ---- Column sanity & basic prep ----
req <- c("settlement_count","mobile_cell_count","stationary_cell_count",
         "motility_prop","day","species")
stopifnot(all(req %in% names(d_spore)))

d_spore <- d_spore |>
  mutate(
    species = as.factor(species),
    total_cell_count = dplyr::coalesce(total_cell_count,
                                       mobile_cell_count + stationary_cell_count),
    total_cell_count = pmax(total_cell_count, 1L),  # avoid log(0)
    day = as.numeric(day)
  )
# If a total column already exists under some name, detect it; otherwise build it
has_total <- any(grepl("^total.*count$", names(d_spore)))
if (!has_total) {
  stopifnot(all(c("mobile_cell_count","stationary_cell_count") %in% names(d_spore)))
  d_spore <- d_spore |>
    mutate(total_cell_count = mobile_cell_count + stationary_cell_count)
} else {
  total_name <- names(d_spore)[grepl("^total.*count$", names(d_spore))][1]
  d_spore <- d_spore |>
    mutate(total_cell_count = !!rlang::sym(total_name))
}

# Final hygiene for the offset and types
d_spore <- d_spore |>
  mutate(
    total_cell_count = pmax(total_cell_count, 1L),
    species = as.factor(species),
    day = as.numeric(day)
  )

# Quick check
d_spore %>% select(mobile_cell_count, stationary_cell_count, total_cell_count) %>% head()

#===============================================================================
# create total_cell_count
library(dplyr)
library(stringr)

# Inspect current names
names(d_spore)


# create total_cell_count (mobile_cell_count + stationary_cell_count)
d_spore$total_cell_count <- pmax(d_spore$mobile_cell_count + d_spore$stationary_cell_count, 1L)


# Final hygiene for the offset and types
d_spore <- d_spore |>
  mutate(
    total_cell_count = pmax(total_cell_count, 1L),
    species = as.factor(species),
    day = as.numeric(day)
  )

# Quick check
d_spore %>% select(mobile_cell_count, stationary_cell_count, total_cell_count) %>% head()

## 1) Map age + make totals
d_spore$day <- as.numeric(d_spore$days_from_start)
d_spore$total_cell_count <- pmax(d_spore$mobile_cell_count + d_spore$stationary_cell_count, 1L)

# 2) Fit NB-GAM with offset(log total) and day smooth
library(mgcv)
m_gam_nb <- gam(
  settlement_count ~ s(motility_prop, k = 5) + s(day, k = 5) + species,
  family = nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML"
)
summary(m_gam_nb)
gam.check(m_gam_nb)

# 3) Age-standardised predictions at day 4 (per 1,000 cells)
mot_vals <- seq(0.30, 0.70, by = 0.05)
tot_med  <- median(d_spore$total_cell_count, na.rm = TRUE)

preds <- do.call(rbind, lapply(levels(as.factor(d_spore$species)), function(sp) {
  nd <- expand.grid(
    motility_prop     = mot_vals,
    day               = 4,
    species           = sp,
    total_cell_count  = tot_med
  )
  p <- predict(m_gam_nb, newdata = nd, type = "link", se.fit = TRUE)
  transform(nd,
            rate_per_1k = exp(p$fit) * 1000,
            lo = exp(p$fit - 1.96*p$se.fit) * 1000,
            hi = exp(p$fit + 1.96*p$se.fit) * 1000)
}))
preds

# refit with larger k and automatic shrinkage

library(mgcv)

m_gam_nb2 <- gam(
  settlement_count ~ 
    s(motility_prop, k = 7) +          # a touch more flexibility
    s(day, k = 10) +                   # give 'day' room to bend
    species,
  family = nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML",
  select = TRUE,    # allows smooths to shrink/penalize extra wiggles
  gamma  = 1.2      # slightly discourages overfitting (optional)
)

summary(m_gam_nb2)
gam.check(m_gam_nb2)


m_gam_nb3 <- mgcv::gam(
  settlement_count ~ 
    s(motility_prop, k = 7) + 
    s(day, k = 15) +             # bump k; unique_days - 1 is a good ceiling
    species,
  family = mgcv::nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML",
  select = TRUE,                 # shrink unneeded wiggles
  gamma  = 1.2
)
gam.check(m_gam_nb3)



m_gam_nb4 <- mgcv::gam(
  settlement_count ~ s(motility_prop, k = 7) + s(day, k = 15) +
    ti(motility_prop, day, k = c(5, 7)),
  family = mgcv::nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML",
  select = TRUE
)
anova(m_gam_nb3, m_gam_nb4, test = "Chisq")


u_day  <- length(unique(d_spore$day))
k_day  <- max(5, min(10, u_day - 1))   # between 5 and 10, but never ≥ unique days
k_day

library(mgcv)

m_gam_nb3 <- gam(
  settlement_count ~ 
    s(motility_prop, k = 7) + 
    s(day, bs = "cr", k = k_day) +   # bs="cr" is stable for 1D; k chosen above
    species,
  family = nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML",
  select = TRUE,
  gamma  = 1.2
)

summary(m_gam_nb3)
gam.check(m_gam_nb3)



# Piecewise linear in day (one kink at day=4)
d_spore$h_day <- pmax(0, d_spore$day - 4)
m_pw <- mgcv::gam(
  settlement_count ~ s(motility_prop, k=7) + day + h_day + species,
  family = nb(link="log"),
  offset = log(total_cell_count),
  data = d_spore, method="REML"
)
anova(m_gam_nb3, m_pw, test="Chisq")   # effects should agree qualitatively

# Or: treat day as a factor (no smooth at all)
m_dayfac <- mgcv::gam(
  settlement_count ~ s(motility_prop, k=7) + factor(day) + species,
  family = nb(link="log"),
  offset = log(total_cell_count),
  data = d_spore, method="REML"
)

#===============================================================================
# Categorical-day model
#===============================================================================
library(mgcv)

# Categorical-day model (no smooth for day)
m_dayfac <- gam(
  settlement_count ~ s(motility_prop, k = 7) + factor(day) + species,
  family = nb(link = "log"),
  offset = log(total_cell_count),
  data   = d_spore,
  method = "REML"
)

summary(m_dayfac)                    # inspect effects
anova(m_gam_nb3, m_dayfac, test="Chisq")  # compare to your smooth-day model

#===============================================================================
# Table of day 4 predictions
#===============================================================================
# Choose which model to use for predictions:
mod <- m_gam_nb3   # (use m_dayfac instead if you want)

mot_vals <- c(0.40, 0.50, 0.60)
tot_med  <- median(d_spore$total_cell_count, na.rm = TRUE)

preds <- do.call(rbind, lapply(levels(as.factor(d_spore$species)), function(sp) {
  nd <- expand.grid(
    motility_prop    = mot_vals,
    day              = 4,               # age-standardised predictions
    species          = sp,
    total_cell_count = tot_med
  )
  p <- predict(mod, newdata = nd, type = "link", se.fit = TRUE)
  transform(nd,
            rate_per_1k = exp(p$fit) * 1000,
            lo          = exp(p$fit - 1.96 * p$se.fit) * 1000,
            hi          = exp(p$fit + 1.96 * p$se.fit) * 1000
  )
}))

# Tidy and save
preds <- preds[order(preds$species, preds$motility_prop), ]
print(preds)
# write.csv(preds, "C:/Motility_Profiling/predictions_day4_per1k_by_species.csv", row.names = FALSE)

#===============================================================================
# Figure to illustrate settlement rate
#===============================================================================

# Choose the model for predictions.
# You found the categorical-day model fits slightly better:
mod <- m_dayfac   # or use m_gam_nb3 if you prefer the smooth-day version

library(dplyr)
library(ggplot2)
library(patchwork)

# Common helpers
mot_seq <- seq(0.30, 0.70, by = 0.01)
tot_med <- median(d_spore$total_cell_count, na.rm = TRUE)
species_lvls <- levels(as.factor(d_spore$species))
day_lvls <- sort(unique(d_spore$day))

## ---------- Panel A: Day-4 curve vs motility (per species) ----------
nd_A <- expand.grid(
  motility_prop    = mot_seq,
  day              = 4,                 # age-standardised at day 4
  species          = species_lvls,
  total_cell_count = tot_med
)

pA_pred <- predict(mod, newdata = nd_A, type = "link", se.fit = TRUE)
nd_A <- nd_A |>
  mutate(rate_per_1k = exp(pA_pred$fit) * 1000,
         lo = exp(pA_pred$fit - 1.96 * pA_pred$se.fit) * 1000,
         hi = exp(pA_pred$fit + 1.96 * pA_pred$se.fit) * 1000)

pA <- ggplot(nd_A, aes(motility_prop, rate_per_1k)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(size = 0.9) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_wrap(~ species, nrow = 1) +
  labs(x = "Motility proportion (day 4)",
       y = "Predicted settlement per 1,000 cells",
       title = "A. Day-4 settlement rate rises steeply with motility") +
  theme_classic(base_size = 11)

## ---------- Panel B: Age curve at motility = 0.5 (per species) ----------
nd_B <- expand.grid(
  motility_prop    = 0.50,              # fixed to show age effect
  day              = day_lvls,
  species          = species_lvls,
  total_cell_count = tot_med
)

pB_pred <- predict(mod, newdata = nd_B, type = "link", se.fit = TRUE)
nd_B <- nd_B |>
  mutate(rate_per_1k = exp(pB_pred$fit) * 1000,
         lo = exp(pB_pred$fit - 1.96 * pB_pred$se.fit) * 1000,
         hi = exp(pB_pred$fit + 1.96 * pB_pred$se.fit) * 1000)

pB <- ggplot(nd_B, aes(day, rate_per_1k, group = species)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(size = 0.9) +
  geom_point(size = 1.6) +
  facet_wrap(~ species, nrow = 1) +
  labs(x = "Days from start",
       y = "Predicted settlement per 1,000 cells",
       title = "B. Even at the same motility (0.5), settlement falls with age") +
  theme_classic(base_size = 11)

## ---------- Combine and save ----------
fig <- pA / pB + plot_layout(heights = c(1, 1.05))

# Sizes: ~170 mm wide x 120 mm tall for journal single-column aesthetics
w_mm <- 170; h_mm <- 120
ggsave("figure_settlement_per1k__motility_day_panels.png", fig, width = w_mm, height = h_mm, units = "mm", dpi = 300)
ggsave("figure_settlement_per1k__motility_day_panels.tiff", fig, width = w_mm, height = h_mm, units = "mm", dpi = 300, compression = "lzw")
ggsave("figure_settlement_per1k__motility_day_panels.pdf",  fig, width = w_mm, height = h_mm, units = "mm")


#===============================================================================
# Developing explanatoiry visual
#===============================================================================

# Use the categorical-day model (fits slightly better)
mod <- m_dayfac   # or m_gam_nb3 if you prefer the smooth-day version

library(dplyr)
library(ggplot2)
library(patchwork)

mot_seq  <- seq(0.30, 0.70, by = 0.01)
tot_med  <- median(d_spore$total_cell_count, na.rm = TRUE)
species_lvls <- levels(as.factor(d_spore$species))

# Helper: make 'day' in newdata match the model's type (factor vs numeric)
coerce_day <- function(x) {
  if (is.factor(d_spore$day)) factor(x, levels = levels(d_spore$day)) else x
}

## ---------- Panel A: Day-4 curve vs motility (per species) ----------
nd_A <- expand.grid(
  motility_prop    = mot_seq,
  day              = coerce_day(4),   # important for m_dayfac
  species          = species_lvls,
  total_cell_count = tot_med
)

pA_pred <- predict(mod, newdata = nd_A, type = "link", se.fit = TRUE)
nd_A <- nd_A |>
  mutate(rate_per_1k = exp(pA_pred$fit) * 1000,
         lo = exp(pA_pred$fit - 1.96 * pA_pred$se.fit) * 1000,
         hi = exp(pA_pred$fit + 1.96 * pA_pred$se.fit) * 1000)

pA <- ggplot(nd_A, aes(motility_prop, rate_per_1k)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +                    # <- use linewidth (fixes warning)
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_wrap(~ species, nrow = 1) +
  labs(x = "Motility proportion (day 4)",
       y = "Predicted settlement per 1,000 cells",
       title = "A. Day-4 settlement rate rises steeply with motility") +
  theme_classic(base_size = 11)

# Show on screen:
print(pA)

## ---------- Panel B: Age curve at motility = 0.5 (per species) ----------
# For m_dayfac, use the observed factor levels; for m_gam_nb3, numeric unique days is fine.
day_vals <- if (is.factor(d_spore$day)) levels(d_spore$day) else sort(unique(d_spore$day))

nd_B <- expand.grid(
  motility_prop    = 0.50,
  day              = coerce_day(day_vals),
  species          = species_lvls,
  total_cell_count = tot_med
)

pB_pred <- predict(mod, newdata = nd_B, type = "link", se.fit = TRUE)
nd_B <- nd_B |>
  mutate(rate_per_1k = exp(pB_pred$fit) * 1000,
         lo = exp(pB_pred$fit - 1.96 * pB_pred$se.fit) * 1000,
         hi = exp(pB_pred$fit + 1.96 * pB_pred$se.fit) * 1000)

pB <- ggplot(nd_B, aes(as.numeric(as.character(day)), rate_per_1k, group = species)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  facet_wrap(~ species, nrow = 1) +
  labs(x = "Days from start",
       y = "Predicted settlement per 1,000 cells",
       title = "B. Even at the same motility (0.5), settlement falls with age") +
  theme_classic(base_size = 11)

# Show on screen:
print(pB)

## ---------- Combine and save (and also show) ----------
fig <- pA / pB + plot_layout(heights = c(1, 1.05))

# Show combined figure on screen:
print(fig)

# Save for the paper
w_mm <- 170; h_mm <- 120
ggsave("C:/Motility_Profiling/plots/figure_settlement_per1k__motility_day_panels.png",
       fig, width = w_mm, height = h_mm, units = "mm", dpi = 300)
ggsave("C:/Motility_Profiling/plots/figure_settlement_per1k__motility_day_panels.tiff",
       fig, width = w_mm, height = h_mm, units = "mm", dpi = 300, compression = "lzw")
ggsave("C:/Motility_Profiling/plots/figure_settlement_per1k__motility_day_panels.pdf",
       fig, width = w_mm, height = h_mm, units = "mm")









#===============================================================================
# 7. Negative-binomial GAM (rate of settlement with offset)
# ==============================================================================
m_gam_nb <- mgcv::gam(
  settlement_count ~
    s(motility_prop, k = 5) +   # allow curvature/threshold
    s(day, k = 5) +             # age control (non-linear)
    species,                    # species shifts
  family = mgcv::nb(link = "log"),
  offset = log(total_cell_count),   # settlement rate per total cells
  data = d_spore,
  method = "REML"
)

summary(m_gam_nb)
mgcv::gam.check(m_gam_nb)   # quick residual/k-index checks



#===============================================================================
# Diagnostics
# DHARMa residual simulation for GAM via predict:
res <- simulateResiduals(m_gam_nb, n = 1000)
plot(res)           # uniformity, dispersion, outliers
testDispersion(res) # over/under-dispersion check


#===============================================================================
# Predictions for the “motility curve” at day 4 (example)
newdat <- expand.grid(
  motility_prop = seq(0.05, 0.85, length.out = 100),
  day = 4,
  species = levels(d_spore$species)[1],        # choose or loop over species
  total_cell_count = median(d_spore$total_cell_count, na.rm = TRUE)
)

pred <- predictions(m_gam_nb, newdata = newdat, type = "link", by = "motility_prop")
# Convert to rate per 1000 cells for interpretability:
pred <- pred %>%
  mutate(rate_per_1k = exp(estimate) * 1000,
         lo = exp(conf.low) * 1000,
         hi = exp(conf.high) * 1000)
head(pred)


#===============================================================================
# Piecewise (“hinge”) model to capture the ~0.5 threshold
library(glmmTMB)
library(splines)

# Hinge at c = 0.5 (change point); adjust if your LOESS suggests ~0.47–0.52
c <- 0.5
d_spore <- d_spore %>%
  mutate(h1 = pmax(0, motility_prop - c))

m_hinge_nb <- glmmTMB(
  settlement_count ~ motility_prop + h1 + day + species +
    offset(log(total_cell_count)),
  family = nbinom2(link = "log"),
  data = d_spore
)

summary(m_hinge_nb)

#===============================================================================
# Diagnostics
simulationOutput <- simulateResiduals(m_hinge_nb, n = 1000)
plot(simulationOutput)
testDispersion(simulationOutput)

#===============================================================================
# Alternative: partial Spearman (settlement ↔ motility, controlling age)
library(ppcor)

# Use raw counts (untransformed) for clarity; or log if you prefer
pc <- pcor.test(
  x = d_spore$motility_prop,
  y = d_spore$settlement_count,
  z = d_spore$day,         # control for age
  method = "spearman"
)
pc # $estimate, $p.value

#===============================================================================
# control both day and species
pc2 <- pcor.test(d_spore$motility_prop,
                 d_spore$settlement_count,
                 d_spore[, c("day","species")], method = "spearman")
pc2



#===============================================================================
# END
# ==============================================================================
