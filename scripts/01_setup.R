# ==============================================================================
# Project:     Kelp Zoospore Motility Profiling
# Script:      01_setup.R
# Author:      Marianne Glascott (with ChatGPT assistant)
# Affiliation: School of Life Sciences, University of Sussex
# Date:        2025-10-31
# Version:     1.0
# ==============================================================================
# Description
#   Create/verify folder tree
#   Load core packages and options
#   Define colour palette, ggplot theme, and scale helpers
#   Define export helpers (PNG/TIFF/PDF)
#   Set up renv (restore/snapshot) and GitHub remote (optional)
#   Provide small utilities (logging, safe package installer)
# Notes
# - Paths are project‑relative using here::here().
# - Figures saved to plots/ (screen/EDA) and figures/ (pub‑quality) in PNG, TIFF, PDF.
# - Git: points origin to https://github.com/MarianneGlascott/motility_profiling.git
# - Designed to work from RStudio Project: C:/Motility_Profiling/Motility_Profiling.Rproj
# ============================================================================
# Contents:
# 1. Repos & lightweight installer
# 2. Anchor project
# 3. Global options
# 4. Packages
# 5. Colour palette & theme
# 6. Figure exporters (PNG, TIFF, PDF + screen)
# 7. renv bootstrap (safe, idempotent)
# 8. Git/GitHub remote (optional; no error if Git absent)
# 9. Logging helper
# ==============================================================================
# 1. Repos & lightweight installer
#===============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

.install_if_missing <- function(pkgs) {
  need <- setdiff(pkgs, rownames(installed.packages()))
  if (length(need)) install.packages(need, quiet = TRUE)
}

.install_if_missing(c("here","fs","glue"))

suppressPackageStartupMessages({
  library(here); library(fs); library(glue)
})
#===============================================================================
# 2. Anchor project (works with or without an .Rproj)
#===============================================================================
if (!fs::file_exists(".here")) fs::file_create(".here")

# ---- 2. Directory constants & create tree -----------------------------------
DIR_ROOT        <- here::here()
DIR_SCRIPTS     <- here::here("scripts")
DIR_DATA_RAW    <- here::here("data_raw")     # already exists per brief
DIR_DATA_OUT    <- here::here("data_processed")
DIR_DOCS        <- here::here("docs")
DIR_FIGURES     <- here::here("figures")
DIR_PLOTS       <- here::here("plots")
DIR_LOGS        <- here::here("logs")
DIR_OUTPUTS     <- here::here("outputs")
DIR_RENV        <- here::here("renv")         # created by renv

# Create any missing project folders
fs::dir_create(c(DIR_SCRIPTS, DIR_DATA_OUT, DIR_DOCS, DIR_FIGURES,
                 DIR_PLOTS, DIR_LOGS, DIR_OUTPUTS))
#===============================================================================
# 3. Global options
#===============================================================================
options(
  scipen = 999,
  digits = 4,
  stringsAsFactors = FALSE,
  readr.show_col_types = FALSE,
  dplyr.summarise.inform = FALSE
)
#===============================================================================
# 4. Packages (respond 'no' to compilation if asked - binary faster!)
#===============================================================================
core_pkgs <- c(
  "tidyverse","janitor","lubridate","scales","ggplot2",
  "ragg","patchwork","ggh4x","ggrepel"
)
.install_if_missing(core_pkgs)

suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(lubridate); library(scales)
  library(ggplot2);   library(ragg);    library(patchwork); library(ggrepel)
})
#===============================================================================
# 5. Colour palette & theme
#===============================================================================
# Kelp palette (named)
kelp_pal <- c(
  kelp_blue  = "#2B6CB0",  # calm ocean blue
  kelp_teal  = "#2C7DA0",  # used previously
  kelp_green = "#2F855A",
  kelp_olive = "#71863F",
  kelp_gold  = "#D4A72C",
  kelp_grey  = "#6E7781",
  kelp_purple= "#7E5AA6"
)

# Discrete scales aligned to class order IMM, NPM, SNM, RPM
scale_fill_kelp  <- function(...) {
  ggplot2::scale_fill_manual(
    values = c(IMM = kelp_pal["kelp_grey"],
               NPM = kelp_pal["kelp_gold"],
               SNM = kelp_pal["kelp_teal"],
               RPM = kelp_pal["kelp_blue"]), ...)
}
scale_colour_kelp <- function(...) {
  ggplot2::scale_colour_manual(
    values = c(IMM = kelp_pal["kelp_grey"],
               NPM = kelp_pal["kelp_gold"],
               SNM = kelp_pal["kelp_teal"],
               RPM = kelp_pal["kelp_blue"]), ...)
}

# Publication‑ready base theme
theme_kelp <- function(base_size = 10, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 6)),
      plot.caption  = element_text(size = rel(0.9), colour = "#555"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

ggplot2::theme_set(theme_kelp())
#===============================================================================
# 6. Figure exporters (PNG, TIFF, PDF + screen)
#===============================================================================
.save_multi <- function(p, dir, base, width, height, units, dpi, bg,
                        formats = c("png", "pdf", "tiff")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  base <- tools::file_path_sans_ext(basename(base))
  has_ragg   <- requireNamespace("ragg", quietly = TRUE)
  has_cairo  <- isTRUE(capabilities("cairo"))
  saved <- character(0)
  for (ext in formats) {
    out <- file.path(dir, paste0(base, ".", ext))
    if (ext == "png" && has_ragg) {
      ggsave(out, p, width = width, height = height, units = units,
             dpi = dpi, bg = bg, device = ragg::agg_png)
    } else if (ext == "tiff" && has_ragg) {
      ggsave(out, p, width = width, height = height, units = units,
             dpi = dpi, bg = bg, device = ragg::agg_tiff)
    } else if (ext == "pdf" && has_cairo) {
      ggsave(out, p, width = width, height = height, units = units,
             dpi = dpi, bg = bg, device = grDevices::cairo_pdf)
    } else if (ext == "tiff") {
      ggsave(out, p, width = width, height = height, units = units,
             dpi = dpi, bg = bg, device = "tiff", compression = "lzw")
    } else {
      ggsave(out, p, width = width, height = height, units = units,
             dpi = dpi, bg = bg)
    }
    message("Saved: ", out); saved <- c(saved, out)
  }
  invisible(saved)
}

export_figure <- function(p, filename, width = 170, height = 120,
                          units = "mm", dpi = 300, bg = "white") {
  .save_multi(p, DIR_FIGURES, filename, width, height, units, dpi, bg)
}

export_plot <- function(p, filename, width = 160, height = 120,
                        units = "mm", dpi = 300, bg = "white") {
  .save_multi(p, DIR_PLOTS, filename, width, height, units, dpi, bg)
}

show_and_save <- function(p, filename, width = 170, height = 120,
                          units = "mm", dpi = 300, bg = "white") {
  if (interactive()) print(p)
  export_plot(p, filename, width = width, height = height,
              units = units, dpi = dpi, bg = bg)
}
#===============================================================================
# 7. renv bootstrap (safe, idempotent)
#===============================================================================
.setup_renv <- function(snapshot = FALSE) {
  .install_if_missing("renv")
  if (!file.exists(here::here("renv.lock"))) {
    renv::init(bare = TRUE)
  } else {
    try(renv::restore(prompt = FALSE), silent = TRUE)
  }
  if (isTRUE(snapshot)) renv::snapshot(prompt = FALSE)
}

# Call now (restore if lock exists)
.setup_renv(snapshot = FALSE)
#===============================================================================
# 8. Git/GitHub remote (optional; no error if Git absent)
#===============================================================================
setup_github_remote <- function(url = "https://github.com/MarianneGlascott/motility_profiling.git") {
  if (nzchar(Sys.which("git"))) {
    .install_if_missing(c("usethis","gitcreds"))
    suppressPackageStartupMessages({ library(usethis) })
    # ensure a repo exists locally
    if (!fs::file_exists(here::here(".git"))) usethis::use_git()
    # set or update origin
    try(usethis::use_git_remote("origin", url, overwrite = TRUE), silent = TRUE)
    # detect branch name
    br <- tryCatch(system("git rev-parse --abbrev-ref HEAD", intern = TRUE), error = function(e) "main")
    if (length(br) == 0 || br %in% c("HEAD","") ) br <- "main"
    # push (won't error fatally if auth missing)
    try(system(glue('git push -u origin {br}'), intern = TRUE), silent = TRUE)
    message("Git remote set to: ", url)
  } else {
    message("Git not found on PATH; skipped GitHub setup.")
  }
}

# Uncomment to run first time
# setup_github_remote()  #first time only

# FIX:
# 1) Make sure Git + packages are ready
install.packages(c("usethis","gert"), quiet = TRUE)

# 2) Initialize a repo if none exists yet
if (!gert::git_is_repo(".")) gert::git_init()

# 3) Ensure at least one commit (needed before pushing)
if (nrow(gert::git_log(max = 1)) == 0) {
  gert::git_add(".")                  # stage everything
  gert::git_commit("Initial commit")  # create first commit
}

# 4) Set the remote and push current branch
url <- "https://github.com/MarianneGlascott/motility_profiling.git"
rem <- gert::git_remote_list()
if (!"origin" %in% rem$name) {
  gert::git_remote_add("origin", url)
} else {
  gert::git_remote_set_url("origin", url)
}

# Detect (or create) current branch and push
br <- gert::git_branch()
current <- br$name[br$head]
if (is.na(current) || current == "") {
  current <- "main"
  if (!"main" %in% br$name) gert::git_branch_create("main")
  gert::git_branch_switch("main")
}

# Push and set upstream (you may get an auth prompt)
try(gert::git_push(remote = "origin", set_upstream = TRUE), silent = TRUE)
message("✅ Remote 'origin' set and push attempted on branch: ", current)

#===============================================================================
# 9. Logging helper
#===============================================================================
log_run <- function(stage, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- glue("[{ts}] {stage}: {msg}\n")
  cat(line, file = file.path(DIR_LOGS, "run.log"), append = TRUE)
  invisible(line)
}

message("✔ 01_setup complete. Root: ", DIR_ROOT)

#===============================================================================
# END
#===============================================================================
