# 03_external_validation.R
# Specificity test in three cohorts independent of TCGA:
#   METABRIC breast cancer (microarray, cBioPortal), CPTAC endometrial carcinoma (RNA-seq, cBioPortal),
#   GSE9891 ovarian carcinoma (microarray, Gene Expression Omnibus).
# Same leave-one-out interferon pathway score, read-outs and controls as in 02_specificity_test_tcga.R;
# genes are matched by symbol. Non-tumor samples and tumors of low malignant potential are excluded.
# Output (output/intermediate): external_validation.csv

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
source(file.path(root, "code", "gene_sets.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(janitor); library(GEOquery) })
set.seed(42); B_BOOT <- 1000

# ---- Genes (symbols; aliases where platforms differ) ------------------------
MARKERS <- list(
  DDX58    = c("DDX58","RIGI"),
  OAS1     = "OAS1",
  OASL     = "OASL",
  GBP2     = "GBP2",
  IRF1     = "IRF1",
  USP18    = "USP18",
  DDX60    = "DDX60",
  LGALS3BP = "LGALS3BP",
  ISG15    = "ISG15",
  IFI6     = "IFI6",
  IFI44    = "IFI44",
  IFI44L   = "IFI44L",
  GBP1     = "GBP1"
)
MECH <- list(NLRC5 = "NLRC5")
IFN <- ifn_symbols("primary")
READ <- c("GZMA","PRF1","CD8A","HLA-A","HLA-B","HLA-C","B2M")
M2_SYM  <- c("CD163","MRC1","MSR1","CD68","CSF1R","VSIG4","MS4A4A","STAB1","MARCO","C1QC")
EXH_SYM <- c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX","ENTPD1","LAYN")
APP_SYM <- c("TAP1","TAP2","TAPBP","PSMB8","PSMB9","PSMB10","ERAP1","ERAP2","CALR","CANX","PDIA3","HLA-A","HLA-B","HLA-C","B2M")
ALLSYM <- unique(c(unlist(MARKERS, use.names=FALSE), unlist(MECH, use.names=FALSE), IFN, READ, M2_SYM, EXH_SYM, APP_SYM))

# ---- Statistics (as in 02_specificity_test_tcga.R) ---------------------------
zsc  <- function(x){ s<-stats::sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s }
getc <- function(df, aliases){ h <- aliases[aliases %in% colnames(df)]
                               if(length(h)) df[[h[1]]] else rep(NA_real_, nrow(df)) }
sp   <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<20) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
part <- function(x,y,Z){ Z<-as.matrix(Z)
  ok<-is.finite(x)&is.finite(y)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  rx<-rank(x[ok]); ry<-rank(y[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  ex<-resid(lm(rx~RZ)); ey<-resid(lm(ry~RZ)); suppressWarnings(cor(ex,ey)) }
incR2 <- function(y,x,Z){ Z<-as.matrix(Z)
  ok<-is.finite(y)&is.finite(x)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  y<-rank(y[ok]); x<-rank(x[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  tryCatch(summary(lm(y~x+RZ))$r.squared - summary(lm(y~RZ))$r.squared, error=function(e) NA_real_) }
# bootstrap 95% CI of the partial correlation (closed form on resampled ranks)
boot_rp <- function(x,y,sc,B=B_BOOT){
  ok<-which(is.finite(x)&is.finite(y)&is.finite(sc)); n<-length(ok); if(n<30) return(c(lo=NA,hi=NA))
  rx<-rank(x[ok]); ry<-rank(y[ok]); rz<-rank(sc[ok]); RP<-numeric(B)
  for(b in 1:B){ s<-sample.int(n,n,replace=TRUE); a<-rx[s]; c2<-ry[s]; d<-rz[s]
    rxy<-cor(a,c2); rxz<-cor(a,d); ryz<-cor(c2,d)
    den<-sqrt((1-rxz^2)*(1-ryz^2)); RP[b]<-if(is.finite(den)&&den>1e-9) (rxy-rxz*ryz)/den else NA_real_ }
  q<-as.numeric(quantile(RP,c(.025,.975),na.rm=TRUE)); c(lo=q[1],hi=q[2]) }

CYT  <- function(df) rowMeans(cbind(zsc(getc(df,"GZMA")), zsc(getc(df,"PRF1"))), na.rm=TRUE)
MHC1 <- function(df) rowMeans(cbind(zsc(getc(df,"HLA-A")), zsc(getc(df,"HLA-B")),
                                    zsc(getc(df,"HLA-C")), zsc(getc(df,"B2M"))), na.rm=TRUE)
# mean z-score of a symbol list without the aliases of the tested gene (leave-one-out)
zmean_sym <- function(df, syms, excl=character(0)){ cols <- setdiff(syms, excl); cols <- cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_, nrow(df))); rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }
M2  <- function(df, excl=character(0)) zmean_sym(df, M2_SYM,  excl)
EXH <- function(df, excl=character(0)) zmean_sym(df, EXH_SYM, excl)
APP <- function(df, excl=character(0)) zmean_sym(df, APP_SYM, excl)
ifn_score <- function(df, excl_aliases=character(0)){
  cols <- setdiff(IFN, excl_aliases); cols <- cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_, nrow(df)))
  rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }

# ---- cBioPortal cohorts (by gene symbol) -----------------------------------
to_wide_sym <- function(md){ md<-janitor::clean_names(as.data.frame(md))
  hs<-grep("hugo|symbol",colnames(md),value=TRUE)[1]
  smp<-grep("^sample_id$|sampleid",colnames(md),value=TRUE)[1]
  val<-grep("^value$",colnames(md),value=TRUE)[1]
  md %>% dplyr::transmute(sample_id=.data[[smp]], sym=.data[[hs]], value=.data[[val]]) %>%
    dplyr::filter(!is.na(sym)) %>% dplyr::distinct(sample_id, sym, .keep_all=TRUE) %>%
    tidyr::pivot_wider(names_from=sym, values_from=value) }

load_cbio <- function(key, study, profile){
  ex <- to_wide_sym(cached_pull(paste0("external_", key, "_", genes_key(ALLSYM)), function(){
    cbio <- cBioPortalData::cBioPortal()
    cBioPortalData::getDataByGenes(cbio, studyId=study, genes=ALLSYM,
                                   by="hugoGeneSymbol", molecularProfileIds=profile) })[[1]])
  # exclude normal tissue samples
  cl <- tryCatch(cached_pull(paste0("external_", key, "_clinical"), function(){
          cbio <- cBioPortalData::cBioPortal()
          as.data.frame(cBioPortalData::clinicalData(cbio, studyId=study)) }), error=function(e) NULL)
  if(!is.null(cl)){
    cl <- janitor::clean_names(cl)
    stc <- grep("sample_type", colnames(cl), value=TRUE)[1]
    sid <- grep("^sample_id$", colnames(cl), value=TRUE)[1]
    if(!is.na(stc) && !is.na(sid)){
      norm <- cl[[sid]][grepl("normal", cl[[stc]], ignore.case=TRUE)]
      if(length(norm)) ex <- ex[!(ex$sample_id %in% norm),]
    }
  }
  ex
}

# ---- GEO cohort --------------------------------------------------------------
load_geo <- function(gse){
  cache_rds <- file.path(path_cache, paste0(tolower(gse), "_expr_", genes_key(ALLSYM), ".rds"))
  if(file.exists(cache_rds)) return(readRDS(cache_rds))
  # series matrix and platform annotation (about 100 MB) are kept in output/cache/geo_raw
  geo_dir <- file.path(path_cache, "geo_raw")
  dir.create(geo_dir, recursive=TRUE, showWarnings=FALSE)
  old_opt <- options(timeout=max(3600, getOption("timeout"))); on.exit(options(old_opt), add=TRUE)
  cat("Downloading:", gse, "from GEO\n")
  g <- tryCatch(GEOquery::getGEO(gse, GSEMatrix=TRUE, getGPL=TRUE, destdir=geo_dir),
                error=function(e){
                  unlink(list.files(geo_dir, full.names=TRUE))
                  stop("GEO download failed (", conditionMessage(e), "). Partial files were removed, please run again.", call.=FALSE) })
  es <- g[[1]]
  X  <- Biobase::exprs(es); fd <- Biobase::fData(es); pd <- Biobase::pData(es)
  symcol <- grep("^Gene.?symbol$|^SYMBOL$", colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  syms <- as.character(fd[[symcol]])
  # probe to gene: per symbol the probe with the highest mean expression ("A /// B" counts for A and B)
  rowmean <- rowMeans(X, na.rm=TRUE)
  sym_list <- lapply(strsplit(syms, "///", fixed=TRUE), trimws)
  long <- data.frame(idx=rep(seq_along(sym_list), lengths(sym_list)),
                     sym=unlist(sym_list), stringsAsFactors=FALSE)
  long <- long[!is.na(long$sym) & long$sym %in% ALLSYM, ]
  pick <- sapply(ALLSYM, function(s){
    idx <- long$idx[long$sym == s]
    if(!length(idx)) return(NA_integer_)
    idx[which.max(rowmean[idx])] })
  found <- names(pick)[!is.na(pick)]
  M <- t(X[pick[!is.na(pick)], , drop=FALSE]); colnames(M) <- found
  E <- as.data.frame(M); E$sample_id <- rownames(M)
  # invasive carcinomas only: exclude tumors of low malignant potential
  txt <- apply(pd, 1, function(r) paste(as.character(r), collapse=" | "))
  lmp <- grepl("LMP|low malignant potential|borderline", txt, ignore.case=TRUE)
  if(any(lmp)) E <- E[!lmp[match(E$sample_id, rownames(pd))], ]
  saveRDS(E, cache_rds); E
}

COHORTS <- list(
  METABRIC   = list(label="METABRIC breast (microarray)",   src="cbio",
                    study="brca_metabric", profile="brca_metabric_mrna"),
  CPTAC_UCEC = list(label="CPTAC endometrial (RNA-seq)",    src="cbio",
                    study="uec_cptac_gdc", profile="uec_cptac_gdc_mrna_seq_tpm"),
  OV_GEO     = list(label="Ovarian GSE9891 (microarray)",   src="geo", gse="GSE9891")
)

# ---- Test per cohort ---------------------------------------------------------
run_cohort <- function(key){
  cfg <- COHORTS[[key]]
  df <- if(cfg$src=="cbio") load_cbio(tolower(key), cfg$study, cfg$profile) else load_geo(cfg$gse)
  cat(key, ": n =", nrow(df), "tumors\n")
  yCYT <- CYT(df); ymhc <- MHC1(df); rows <- list(); k <- 1
  # genes of each read-out, removed from the pathway score when testing against that read-out
  RO_SYM <- list(CYT=c("GZMA","PRF1"), M2=M2_SYM, EXH=EXH_SYM, APP=APP_SYM, MHC1=c("HLA-A","HLA-B","HLA-C","B2M"))
  add_row <- function(g, rn, typ, x, yy, sc){ ci <- boot_rp(x,yy,sc)
    data.frame(cohort=key, label=cfg$label, gene=g, readout=rn, type=typ,
      n=sum(is.finite(x)&is.finite(yy)&is.finite(sc)),
      rho_raw=round(sp(x,yy),3), rho_partial=round(part(x,yy,sc),3),
      rp_lo=round(ci["lo"],3), rp_hi=round(ci["hi"],3),
      unique_R2=round(incR2(yy,x,sc),3), row.names=NULL) }
  for(g in names(MARKERS)){
    al <- MARKERS[[g]]; x <- getc(df, al)
    ro_vecs <- list(CYT=yCYT, M2=M2(df, al), EXH=EXH(df, al), APP=APP(df, al))
    for(rn in names(ro_vecs)){
      sc <- ifn_score(df, excl_aliases=c(al, RO_SYM[[rn]]))
      rows[[k]] <- add_row(g, rn, "candidate", x, ro_vecs[[rn]], sc); k<-k+1 }
  }
  # controls against both main read-outs, NLRC5 also against MHC class I
  for(rn in c("CYT","EXH")){ yy <- if(rn=="CYT") yCYT else EXH(df)
    rows[[k]] <- add_row("CD8A",  rn, "positive control",  getc(df,"CD8A"),    yy, ifn_score(df, c("CD8A",  RO_SYM[[rn]]))); k<-k+1
    rows[[k]] <- add_row("NLRC5", rn, "mechanism control", getc(df,MECH$NLRC5), yy, ifn_score(df, c(MECH$NLRC5, RO_SYM[[rn]]))); k<-k+1 }
  rows[[k]] <- add_row("NLRC5", "MHC1", "mechanism control", getc(df,MECH$NLRC5), ymhc, ifn_score(df, c(MECH$NLRC5, RO_SYM$MHC1)))
  dplyr::bind_rows(rows)
}

V <- dplyr::bind_rows(lapply(names(COHORTS), run_cohort))
write.csv(V, file.path(path_inter,"external_validation.csv"), row.names=FALSE)
cat("Done: 03_external_validation.R\n")
