# =============================================================================
# run_all.R  Runs the complete two-stage analysis
#
# Usage: open Revision_TwoStage_20261001.Rproj (or setwd() to that folder) and
#        source("R/run_all.R", encoding = "UTF-8")
#        or, from a terminal in that folder:  Rscript R/run_all.R
# Packages: cluster, ggplot2, dplyr, tidyr, readr, jsonlite, ordinal
# Default sizes: 1000 coefficient bootstraps, 100 splits, 1000 random sets (environment
# variables PM_BOOT, PM_SPLIT, PM_RANDOM override them). Running time is about 3 minutes. Results go to output/tables, output/figures, output/objects.
# =============================================================================
if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args_all <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", args_all, value = TRUE)
ROOT <- if (length(f)) normalizePath(file.path(dirname(sub("^--file=", "", f[1])), ".."), winslash = "/", mustWork = TRUE) else
        normalizePath(if (file.exists("data/phenotypes.csv")) "." else "..", winslash = "/", mustWork = TRUE)
started <- Sys.time()
for (script in c("00_config.R", "01_prepare.R", "02_reproducibility.R", "03_structure.R", "04_validation.R",
                 "05_priority.R", "06_figures.R", "07_export.R")) {
  message("\n== ", script, "  (", format(Sys.time(), "%H:%M:%S"), ")")
  source(file.path(ROOT, "R", script), encoding = "UTF-8")
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
saveRDS(list(data = dat, analysis = main, clustering = base_run, runs_summary = sens, core = core_ids, shinycore = sc), file.path(OUT, "objects", "analysis_objects.rds"))
message("\nCompleted in ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " minutes")
