# Installs the R packages required by the analysis (CRAN and Bioconductor), if missing.
# Versions used for the paper are listed in sessionInfo.txt.

cran <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "janitor", "readxl", "msigdbr", "BiocManager")
miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)

bioc <- c("cBioPortalData", "GSVA", "GEOquery", "Biobase")
miss <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) BiocManager::install(miss, ask = FALSE, update = FALSE)
