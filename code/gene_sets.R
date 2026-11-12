# Interferon gene sets from MSigDB (msigdbr), sourced after setup.R.
#   primary        REACTOME_INTERFERON_SIGNALING (interferon pathway score)
#   hallmark_ifng  HALLMARK_INTERFERON_GAMMA_RESPONSE (sensitivity analysis)
#   hallmark_ifna  HALLMARK_INTERFERON_ALPHA_RESPONSE (sensitivity analysis)
# Provides ifn_symbols(key), ifn_entrez(key), the table IFN_GENESETS and genes_key().

IFN_SETS <- c(primary       = "REACTOME_INTERFERON_SIGNALING",
              hallmark_ifng = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              hallmark_ifna = "HALLMARK_INTERFERON_ALPHA_RESPONSE")

MSIGDBR_VERSION_PAPER <- "26.1.1"
if (as.character(utils::packageVersion("msigdbr")) != MSIGDBR_VERSION_PAPER)
  warning("msigdbr ", as.character(utils::packageVersion("msigdbr")), " is installed, the paper used ",
          MSIGDBR_VERSION_PAPER, ". Other versions can contain different gene sets.", call. = FALSE)

# Eight-character checksum of a gene list, used in cache file names
genes_key <- function(x){
  tf <- tempfile(); writeLines(sort(unique(as.character(x))), tf)
  k  <- substr(unname(tools::md5sum(tf)), 1, 8); unlink(tf); k }

.msig_std <- function(d){ d <- as.data.frame(d)
  gs <- grep("^gs_name$", colnames(d), value=TRUE)[1]
  sy <- grep("^gene_symbol$", colnames(d), value=TRUE)[1]
  en <- grep("^ncbi_gene$|^entrez_gene$", colnames(d), value=TRUE)[1]
  if(is.na(en)) en <- grep("ncbi_gene|entrez", colnames(d), value=TRUE)[1]
  data.frame(gs_name=as.character(d[[gs]]), symbol=as.character(d[[sy]]),
             entrez=suppressWarnings(as.integer(d[[en]])), stringsAsFactors=FALSE) }
# msigdbr >= 10 uses the argument "collection", older versions "category"
.msig_get <- function(coll) .msig_std(tryCatch(msigdbr::msigdbr(species="Homo sapiens", collection=coll),
                                               error=function(e) msigdbr::msigdbr(species="Homo sapiens", category=coll)))

.ifn_cache <- file.path(path_cache, paste0("msigdb_ifn_sets_", as.character(utils::packageVersion("msigdbr")), ".rds"))
if(file.exists(.ifn_cache)){
  .IFNSETS <- readRDS(.ifn_cache)
} else {
  .IFNSETS <- rbind(.msig_get("H"), .msig_get("C2"))
  .IFNSETS <- .IFNSETS[.IFNSETS$gs_name %in% IFN_SETS, ]
  saveRDS(.IFNSETS, .ifn_cache)
}
ifn_symbols <- function(key="primary"){ s <- unique(.IFNSETS$symbol[.IFNSETS$gs_name == IFN_SETS[[key]]]); s[!is.na(s) & nzchar(s)] }
ifn_entrez  <- function(key="primary"){ e <- unique(.IFNSETS$entrez[.IFNSETS$gs_name == IFN_SETS[[key]]]); e[!is.na(e)] }
for(.k in names(IFN_SETS)) if(length(ifn_symbols(.k)) < 50)
  stop("Gene set ", IFN_SETS[[.k]], " was not loaded completely. Check msigdbr and delete ", .ifn_cache, call.=FALSE)
IFN_GENESETS <- do.call(rbind, lapply(names(IFN_SETS), function(k){
  d <- unique(.IFNSETS[.IFNSETS$gs_name == IFN_SETS[[k]], c("symbol","entrez")])
  data.frame(gene_set=k, gs_name=IFN_SETS[[k]], d[order(d$symbol), ], row.names=NULL, stringsAsFactors=FALSE) }))
