# Runs the complete analysis in order, each script in a fresh R session.
# From the repository root:  Rscript run_all.R
# Writes a run log (output/run_log.txt) and the session information (output/sessionInfo.txt).

root <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(root, "code", "setup.R"))) stop("Run run_all.R from the repository root.", call. = FALSE)
Sys.setenv(IFN_SPEC_ROOT = root)
source(file.path(root, "requirements.R"))

scripts <- c("01_threshold_and_score_construction.R",
             "02_specificity_test_tcga.R",
             "03_external_validation.R",
             "04_calibration_simulation.R",
             "05_proliferation_control.R",
             "06_figures_and_tables.R")
START_AT <- 1   # to resume after an interruption, e.g. 3 starts with 03_external_validation.R

rscript  <- file.path(R.home("bin"), "Rscript")
log_file <- file.path(root, "output", "run_log.txt")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logline <- function(...) { l <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...); cat(l, "\n"); cat(l, "\n", file = log_file, append = TRUE) }

if (START_AT == 1) cat("", file = log_file)
logline("Start, ", R.version.string)
t_all <- Sys.time()
for (s in scripts[START_AT:length(scripts)]) {
  cat("\n==========", s, "==========\n")
  t0 <- Sys.time()
  status <- system2(rscript, shQuote(file.path(root, "code", s)))
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (status != 0) { logline("ERROR  ", s, "  (", mins, " min)"); stop("Stopped in ", s, ".", call. = FALSE) }
  logline("OK     ", s, "  (", mins, " min)")
}
logline("Finished, total ", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min")

pkgs <- c("cBioPortalData", "GSVA", "GEOquery", "Biobase", "msigdbr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "patchwork", "janitor", "readxl")
invisible(lapply(pkgs, requireNamespace, quietly = TRUE))
writeLines(capture.output(sessionInfo()), file.path(root, "output", "sessionInfo.txt"))
