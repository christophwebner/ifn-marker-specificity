# 01_threshold_and_score_construction.R
# (A) Raw and partial correlation with cytolytic activity (CYT) for every gene of the Reactome
#     interferon signaling set, its expression-matched non-interferon partner and the thirteen
#     published markers; the 95th percentile of |partial| of the matched genes is the empirical threshold.
# (B) Raw and partial correlation of the markers and controls under six constructions of the
#     interferon pathway score (CYT, MHC class I, antigen processing and presentation).
# The pathway score never contains the tested gene (leave-one-out) or the genes of the read-out.
# Output (output/intermediate): ifn_gene_sets.csv, background_correlations.csv,
#   empirical_threshold.csv, score_construction.csv

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
source(file.path(root, "code", "gene_sets.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(janitor) })
set.seed(42)

# ---- Helpers ---------------------------------------------------------------
zsc <- function(x){ s<-stats::sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s }
# partial Spearman correlation: Pearson correlation of rank residuals after regression on the ranked covariate
part <- function(x,y,z){ ok<-is.finite(x)&is.finite(y)&is.finite(z); if(sum(ok)<20) return(NA_real_)
  rx<-rank(x[ok]); ry<-rank(y[ok]); rz<-rank(z[ok]); ex<-resid(lm(rx~rz)); ey<-resid(lm(ry~rz)); suppressWarnings(cor(ex,ey)) }
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<20) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
to_wide_sym <- function(md){ md<-janitor::clean_names(as.data.frame(md))
  hs<-grep("hugo|symbol",colnames(md),value=TRUE)[1]; smp<-grep("^sample_id$|sampleid",colnames(md),value=TRUE)[1]; val<-grep("^value$",colnames(md),value=TRUE)[1]
  md %>% dplyr::transmute(sample_id=.data[[smp]], sym=.data[[hs]], value=.data[[val]]) %>%
    dplyr::filter(!is.na(sym)) %>% dplyr::distinct(sample_id,sym,.keep_all=TRUE) %>% tidyr::pivot_wider(names_from=sym,values_from=value) }
getcol <- function(df,cands){ h<-cands[cands %in% colnames(df)]; if(length(h)) df[[h[1]]] else rep(NA_real_,nrow(df)) }

# ---- Gene sets -------------------------------------------------------------
write.csv(IFN_GENESETS, file.path(path_inter,"ifn_gene_sets.csv"), row.names=FALSE)
# All analyses in this script use CYT (GZMA, PRF1); these genes are removed from the interferon sets.
CYT_SYMS <- c("GZMA","PRF1")
ifn_p <- setdiff(ifn_symbols("primary"),       CYT_SYMS)
ifn_g <- setdiff(ifn_symbols("hallmark_ifng"), CYT_SYMS)
ifn_a <- setdiff(ifn_symbols("hallmark_ifna"), CYT_SYMS)

# ---- Frozen design: background gene pool and expression-matched pairs ------
# The pool of 2000 non-interferon genes and the matched pairs are fixed inputs (data/).
# REBUILD_DESIGN <- TRUE draws a new pool and rebuilds the pairs with the same procedure,
# which changes the background and the threshold.
REBUILD_DESIGN <- FALSE
bg_file    <- file.path(path_data, "background_gene_pool.csv")
pairs_file <- file.path(path_data, "matched_gene_pairs.csv")
if(!REBUILD_DESIGN){
  bg_pool <- read.csv(bg_file, stringsAsFactors=FALSE)$gene_symbol
  bg_pool <- bg_pool[!is.na(bg_pool) & nzchar(bg_pool)]
} else {
  allmh <- tryCatch(as.data.frame(msigdbr::msigdbr(species="Homo sapiens")), error=function(e) NULL)
  universe <- if(!is.null(allmh)) unique(allmh[[grep("gene_symbol",colnames(allmh),value=TRUE)[1]]]) else character(0)
  universe <- setdiff(universe, c(ifn_p, ifn_g, ifn_a)); universe <- universe[!is.na(universe)&nzchar(universe)]
  bg_pool  <- sample(universe, min(2000, length(universe)))
  write.csv(data.frame(gene_symbol=bg_pool), bg_file, row.names=FALSE)
}

# Thirteen published markers (as in 02_specificity_test_tcga.R)
marker_alias <- list("RIG-I"=c("DDX58","RIGI"),"OAS1"="OAS1","OASL"="OASL","GBP2"="GBP2",
                     "IRF1"="IRF1","USP18"="USP18","DDX60"="DDX60","LGALS3BP"="LGALS3BP",
                     "GBP1"="GBP1","ISG15"="ISG15","IFI6"="IFI6","IFI44"="IFI44","IFI44L"="IFI44L")
extra  <- unique(c(unlist(marker_alias,use.names=FALSE),"CD8A","GZMA","PRF1"))
RO_SETS <- list(
  CYT  = c("GZMA","PRF1"),
  MHC1 = c("HLA-A","HLA-B","HLA-C","B2M"),
  APP  = c("TAP1","TAP2","TAPBP","PSMB8","PSMB9","PSMB10","ERAP1","ERAP2","CALR","CANX","PDIA3",
           "HLA-A","HLA-B","HLA-C","B2M"))
allsym <- unique(c(ifn_p, ifn_g, ifn_a, bg_pool, extra, "NLRC5", unlist(RO_SETS, use.names=FALSE)))

# ---- TCGA PanCancer Atlas expression (cBioPortal, RSEM) --------------------
STUDIES <- c(UCEC="ucec_tcga_pan_can_atlas_2018", OV="ov_tcga_pan_can_atlas_2018", CESC="cesc_tcga_pan_can_atlas_2018", BRCA="brca_tcga_pan_can_atlas_2018")
pull_sym <- function(nm, syms){ study<-STUDIES[nm]; prof<-paste0(study,"_rna_seq_v2_mrna")
  to_wide_sym(cached_pull(paste0("tcga_",tolower(nm),"_background_",genes_key(syms)), function(){ cbio<-cBioPortalData::cBioPortal()
    cBioPortalData::getDataByGenes(cbio, studyId=study, genes=unique(syms), by="hugoGeneSymbol", molecularProfileIds=prof) })[[1]]) }
COH <- lapply(names(STUDIES), function(nm) pull_sym(nm, allsym)); names(COH)<-names(STUDIES)

# Greedy matching without replacement on mean log2 expression in UCEC, most highly expressed genes first
uc <- COH[["UCEC"]]
if(!REBUILD_DESIGN){
  pairs <- read.csv(pairs_file, stringsAsFactors=FALSE)
} else {
prog_present <- intersect(ifn_p, colnames(uc))
bg_present   <- setdiff(intersect(bg_pool, colnames(uc)), c(ifn_p, ifn_g, ifn_a, extra))
memean <- function(df,g) mean(df[[g]], na.rm=TRUE)
lg <- function(x) log2(pmax(x, 0) + 1)
bg_mean_raw <- sapply(bg_present,  function(g) memean(uc, g))
bg_meanL    <- lg(bg_mean_raw)
prog_meanL  <- sapply(prog_present, function(g) lg(memean(uc, g)))
ord <- names(sort(prog_meanL, decreasing=TRUE))
pairs <- data.frame(ifn=character(), ctrl=character(), ifn_mean=numeric(), ctrl_mean=numeric(),
                    ifn_log2=numeric(), ctrl_log2=numeric(), stringsAsFactors=FALSE)
used <- character(0)
for(g in ord){
  mg <- prog_meanL[[g]]
  avail <- bg_present[!(bg_present %in% used)]
  if(!length(avail)) break
  cand <- avail[ which.min(abs(bg_meanL[avail] - mg)) ]
  pairs <- rbind(pairs, data.frame(ifn=g, ctrl=cand, ifn_mean=memean(uc, g), ctrl_mean=unname(bg_mean_raw[[cand]]),
                                   ifn_log2=mg, ctrl_log2=unname(bg_meanL[[cand]]), stringsAsFactors=FALSE))
  used <- c(used, cand)
}
pairs$abs_diff_log2 <- abs(pairs$ifn_log2 - pairs$ctrl_log2)
write.csv(pairs[order(-pairs$ifn_mean), ], pairs_file, row.names=FALSE)
}

# ---- Scores ----------------------------------------------------------------
set_score <- function(df, set, excl=character(0)){ cols <- intersect(setdiff(set, excl), colnames(df))
  if(!length(cols)) return(rep(NA_real_,nrow(df))); rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }
pca1_score <- function(df, set){ cols<-intersect(set, colnames(df)); Z<-sapply(cols, function(c) zsc(df[[c]])); Z[is.na(Z)]<-0
  pc<-tryCatch(prcomp(Z, center=FALSE, scale.=FALSE), error=function(e) NULL); if(is.null(pc)) return(rep(NA_real_,nrow(df)))
  s<-pc$x[,1]; if(suppressWarnings(cor(s, rowMeans(Z)))<0) s<--s; s }
cyt_of <- function(df) rowMeans(cbind(zsc(getcol(df,"GZMA")), zsc(getcol(df,"PRF1"))), na.rm=TRUE)

# ---- (A) Pathway genes, matched genes and markers vs. CYT ------------------
rows<-list(); k<-1
for(nm in names(COH)){ df<-COH[[nm]]; CYTv<-cyt_of(df)
  scF<-set_score(df, ifn_p)                                    # full score for non-member genes
  for(g in intersect(ifn_p, colnames(df))){ sc<-set_score(df, ifn_p, excl=g); rr<-sp(df[[g]],CYTv); rp<-part(df[[g]],CYTv,sc)
    rows[[k]]<-data.frame(cohort=nm, gene=g, role="IFN-programme gene", raw=rr, partial=rp, reduction=rr-rp); k<-k+1 }
  for(c in intersect(pairs$ctrl, colnames(df))){ rr<-sp(df[[c]],CYTv); rp<-part(df[[c]],CYTv,scF)
    rows[[k]]<-data.frame(cohort=nm, gene=c, role="matched non-IFN gene", raw=rr, partial=rp, reduction=rr-rp); k<-k+1 }
  for(mn in names(marker_alias)){ mv<-getcol(df, marker_alias[[mn]]); if(all(is.na(mv))) next
    sc<-set_score(df, ifn_p, excl=marker_alias[[mn]]); rr<-sp(mv,CYTv); rp<-part(mv,CYTv,sc)
    rows[[k]]<-data.frame(cohort=nm, gene=mn, role="published marker", raw=rr, partial=rp, reduction=rr-rp); k<-k+1 } }
R <- dplyr::bind_rows(rows)
# median over the four cohorts per gene
G <- R %>% dplyr::group_by(gene, role) %>% dplyr::summarise(raw=median(raw,na.rm=TRUE), partial=median(partial,na.rm=TRUE),
        reduction=median(reduction,na.rm=TRUE), .groups="drop")
write.csv(G, file.path(path_inter,"background_correlations.csv"), row.names=FALSE)

# Empirical threshold: 95th percentile of |partial| of the expression-matched non-interferon genes
THRESHOLD_PCT <- 0.95
ctrlP_abs <- abs(G$partial[G$role=="matched non-IFN gene"])
threshold <- as.numeric(quantile(ctrlP_abs, THRESHOLD_PCT, na.rm=TRUE))
write.csv(data.frame(percentile=THRESHOLD_PCT, threshold=round(threshold,4), n_matched=sum(!is.na(ctrlP_abs))),
          file.path(path_inter,"empirical_threshold.csv"), row.names=FALSE)
cat(sprintf("Empirical threshold |partial| = %.4f (n = %d matched genes)\n", threshold, sum(!is.na(ctrlP_abs))))

# ---- (B) Six constructions of the interferon pathway score -----------------
if(!requireNamespace("GSVA", quietly=TRUE)) stop("Package GSVA is required (see requirements.R).", call.=FALSE)
# Expression matrix (genes x samples) of all downloaded genes for GSVA and ssGSEA; missing values set to the gene median
expr_matrix <- function(df){
  keep <- vapply(df, is.numeric, logical(1)); M <- t(as.matrix(df[, keep, drop=FALSE])); storage.mode(M)<-"double"
  for(i in seq_len(nrow(M))){ r<-M[i,]; if(anyNA(r)){ m<-stats::median(r,na.rm=TRUE); if(is.na(m)) m<-0; r[is.na(r)]<-m; M[i,]<-r } }
  M }
# GSVA or ssGSEA scores for several gene sets in one call (samples x gene sets); identical sets are scored once
gsva_scores <- function(M, method, sets){
  sets <- lapply(sets, function(s) intersect(s, rownames(M)))
  out  <- matrix(NA_real_, ncol(M), length(sets), dimnames=list(NULL, names(sets)))
  ok   <- vapply(sets, length, integer(1)) >= 5
  key  <- vapply(sets, function(s) paste(sort(s), collapse="|"), character(1))
  uniq <- ok & !duplicated(key)
  par  <- if(method=="gsva") GSVA::gsvaParam(M, sets[uniq]) else GSVA::ssgseaParam(M, sets[uniq], normalize=FALSE)
  sc   <- GSVA::gsva(par, verbose=FALSE)
  for(s in intersect(rownames(sc), names(sets)[uniq])) out[, s] <- as.numeric(sc[s, ])
  for(s in names(sets)[ok & !uniq]) out[, s] <- out[, names(sets)[match(key[[s]], key)]]
  out }

# Tested genes per read-out: the thirteen markers plus the matching control (CYT: CD8A; MHC class I and APP: NLRC5)
CTRL_BY_RO <- list(CYT=list("CD8A"="CD8A"), MHC1=list("NLRC5"="NLRC5"), APP=list("NLRC5"="NLRC5"))
SETS0 <- list(primary=ifn_symbols("primary"), hallmark_ifng=ifn_symbols("hallmark_ifng"), hallmark_ifna=ifn_symbols("hallmark_ifna"))
ro_score <- function(df, r) rowMeans(sapply(RO_SETS[[r]], function(g) zsc(getcol(df, g))), na.rm=TRUE)
robrows<-list(); j<-1
for(nm in names(COH)){ df<-COH[[nm]]; M<-expr_matrix(df)
  jobs <- do.call(rbind, lapply(names(RO_SETS), function(r) data.frame(readout=r,
            gene=c(names(marker_alias), names(CTRL_BY_RO[[r]])), stringsAsFactors=FALSE)))
  jobs$key <- paste(jobs$readout, jobs$gene, sep="|")
  alias_of <- function(g) c(marker_alias, CTRL_BY_RO$CYT, CTRL_BY_RO$MHC1)[[g]]
  loo_sets <- setNames(lapply(seq_len(nrow(jobs)), function(i)
                 setdiff(SETS0$primary, c(alias_of(jobs$gene[i]), RO_SETS[[jobs$readout[i]]]))), jobs$key)
  GS <- gsva_scores(M, "gsva", loo_sets); SS <- gsva_scores(M, "ssgsea", loo_sets)
  for(i in seq_len(nrow(jobs))){ r<-jobs$readout[i]; mn<-jobs$gene[i]; al<-alias_of(mn)
    mv<-getcol(df, al); if(all(is.na(mv))) next
    Y <- ro_score(df, r); ex <- c(al, RO_SETS[[r]])
    robrows[[j]]<-data.frame(cohort=nm, readout=r, gene=mn, role=if(mn %in% names(marker_alias)) "marker" else "control",
      raw                 = sp(mv,Y),
      primary_meanz       = part(mv,Y,set_score(df, SETS0$primary, excl=ex)),
      primary_pc1         = part(mv,Y,pca1_score(df, setdiff(SETS0$primary, ex))),
      primary_gsva        = part(mv,Y,GS[, jobs$key[i]]),
      primary_ssgsea      = part(mv,Y,SS[, jobs$key[i]]),
      hallmark_ifng_meanz = part(mv,Y,set_score(df, SETS0$hallmark_ifng, excl=ex)),
      hallmark_ifna_meanz = part(mv,Y,set_score(df, SETS0$hallmark_ifna, excl=ex))); j<-j+1 } }
ROB_COH <- dplyr::bind_rows(robrows)
VAR_COLS <- c("raw","primary_meanz","primary_pc1","primary_gsva","primary_ssgsea","hallmark_ifng_meanz","hallmark_ifna_meanz")
ROB <- ROB_COH %>% dplyr::group_by(readout, gene, role) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(VAR_COLS), ~ round(median(.x, na.rm=TRUE), 3)), .groups="drop")
write.csv(ROB, file.path(path_inter,"score_construction.csv"), row.names=FALSE)
cat("Done: 01_threshold_and_score_construction.R\n")
