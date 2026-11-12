# Shared setup, sourced at the start of every analysis script:
# repository paths and a cache for data downloads.

root <- Sys.getenv("IFN_SPEC_ROOT", unset = "")
if (root == "") root <- normalizePath(".")
path_data   <- file.path(root, "data")
path_output <- file.path(root, "output")
path_inter  <- file.path(path_output, "intermediate")
path_cache  <- file.path(path_output, "cache")
if (!file.exists(file.path(path_data, "background_gene_pool.csv")))
  stop("Repository root not found. Run the scripts from the repository root or set IFN_SPEC_ROOT.", call. = FALSE)
for (d in c(path_output, path_inter, path_cache)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Returns the cached result if it exists, otherwise runs pull_function() once and caches its result.
cached_pull <- function(cache_name, pull_function) {
  cache_file <- file.path(path_cache, paste0(cache_name, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  cat("Downloading:", cache_name, "\n")
  result <- pull_function()
  saveRDS(result, cache_file)
  result
}
