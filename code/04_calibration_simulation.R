# 04_calibration_simulation.R
# Semi-synthetic calibration of the specificity test in the four TCGA cohorts.
# Synthetic genes G = a*S + b*U + E are embedded in the observed data, with S the interferon
# pathway score, U the part of cytolytic activity (CYT) independent of S and E noise of fixed
# variance; a sets the co-regulation with the pathway, b the gene-specific association with CYT.
# Each synthetic gene is tested with the same rank-based partial correlation as the real genes.
# The real markers and controls are placed on the same axes.
# Output (output/intermediate): simulation_grid.csv, simulation_real_genes.csv

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
source(file.path(root, "code", "gene_sets.R"))
suppressPackageStartupMessages(library(dplyr))
set.seed(42)

# ---- Helpers ---------------------------------------------------------------
zsc <- function(x){ s<-stats::sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s }
part <- function(x,y,z){ ok<-is.finite(x)&is.finite(y)&is.finite(z); if(sum(ok)<20) return(NA_real_)
  rx<-rank(x[ok]); ry<-rank(y[ok]); rz<-rank(z[ok]); ex<-resid(lm(rx~rz)); ey<-resid(lm(ry~rz)); suppressWarnings(cor(ex,ey)) }
to_wide_sym <- function(md){ md<-janitor::clean_names(as.data.frame(md))
  hs<-grep("hugo|symbol",colnames(md),value=TRUE)[1]; smp<-grep("^sample_id$|sampleid",colnames(md),value=TRUE)[1]; val<-grep("^value$",colnames(md),value=TRUE)[1]
  md %>% dplyr::transmute(sample_id=.data[[smp]], sym=.data[[hs]], value=.data[[val]]) %>%
    dplyr::filter(!is.na(sym)) %>% dplyr::distinct(sample_id,sym,.keep_all=TRUE) %>% tidyr::pivot_wider(names_from=sym,values_from=value) }
getcol <- function(df,cands){ h<-cands[cands %in% colnames(df)]; if(length(h)) df[[h[1]]] else rep(NA_real_,nrow(df)) }

IFN_SYM <- ifn_symbols("primary")
ifn_score <- function(df, excl=character(0)){ cols<-setdiff(intersect(IFN_SYM, colnames(df)), excl); if(!length(cols)) return(rep(NA_real_,nrow(df)))
  rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }

RETAINED <- c("IRF1","GBP1","GBP2")
PROXY    <- c("DDX58","OAS1","OASL","USP18","DDX60","LGALS3BP","ISG15","IFI6","IFI44","IFI44L")
CTRL     <- c("CD8A","NLRC5")
MARKER_SYMS <- c(RETAINED, PROXY, CTRL)

# ---- TCGA PanCancer Atlas expression (cBioPortal, RSEM) --------------------
STUDIES <- c(UCEC="ucec_tcga_pan_can_atlas_2018", OV="ov_tcga_pan_can_atlas_2018",
             CESC="cesc_tcga_pan_can_atlas_2018", BRCA="brca_tcga_pan_can_atlas_2018")
need_sym <- unique(c(IFN_SYM, "GZMA","PRF1", MARKER_SYMS, "RIGI"))   # RIGI: current HGNC symbol of DDX58
pull_sym <- function(nm){ study<-STUDIES[nm]; prof<-paste0(study,"_rna_seq_v2_mrna")
  to_wide_sym(cached_pull(paste0("tcga_",tolower(nm),"_simulation_",genes_key(need_sym)), function(){ cbio<-cBioPortalData::cBioPortal()
    cBioPortalData::getDataByGenes(cbio, studyId=study, genes=need_sym, by="hugoGeneSymbol", molecularProfileIds=prof) })[[1]]) }
COH <- lapply(names(STUDIES), function(nm) pull_sym(nm)); names(COH)<-names(STUDIES)

THRESHOLD <- as.numeric(read.csv(file.path(path_inter,"empirical_threshold.csv"))$threshold[1])

# ---- Simulation --------------------------------------------------------------
A_GRID <- c(0, 0.3, 0.6, 0.9, 1.2)      # co-regulation with the pathway
B_GRID <- c(0, 0.15, 0.3, 0.5, 0.8)     # gene-specific association with CYT
NREP   <- 200                           # synthetic genes per (a, b) and cohort

CYT_SYMS <- c("GZMA","PRF1")   # read-out genes are never part of the pathway score
sim_cohort <- function(df){
  S <- ifn_score(df, excl=CYT_SYMS)
  Y <- rowMeans(cbind(zsc(getcol(df,"GZMA")), zsc(getcol(df,"PRF1"))), na.rm=TRUE)
  ok <- is.finite(S) & is.finite(Y); S<-S[ok]; Y<-Y[ok]; n<-length(S)
  zS <- zsc(S); zY <- zsc(Y)
  U  <- resid(lm(zY ~ zS)); zU <- zsc(U)
  out <- list(); j <- 1
  for(a in A_GRID) for(b in B_GRID){
    parts <- numeric(NREP); coregs <- numeric(NREP)
    for(r in seq_len(NREP)){
      E <- rnorm(n)
      G <- a*zS + b*zU + E
      parts[r]  <- part(G, Y, S)
      coregs[r] <- abs(suppressWarnings(cor(G, S, method="spearman")))   # realized co-regulation
    }
    out[[j]] <- data.frame(a=a, b=b,
      median_partial = median(parts, na.rm=TRUE),
      coreg = median(coregs, na.rm=TRUE),
      pass_rate = mean(abs(parts) >= THRESHOLD, na.rm=TRUE),
      n_rep = sum(is.finite(parts)))
    j <- j+1
  }
  dplyr::bind_rows(out)
}
SIM <- dplyr::bind_rows(lapply(names(COH), function(nm){ d<-sim_cohort(COH[[nm]]); d$cohort<-nm; d }))
write.csv(SIM, file.path(path_inter,"simulation_grid.csv"), row.names=FALSE)

# ---- Real genes: co-regulation with the pathway score and partial correlation with CYT --
ALIASES <- list(DDX58=c("DDX58","RIGI"))
real_one <- function(sym){
  xs<-c(); ys<-c(); cands <- if(sym %in% names(ALIASES)) ALIASES[[sym]] else sym
  for(nm in names(COH)){ df<-COH[[nm]]
    Y <- rowMeans(cbind(zsc(getcol(df,"GZMA")), zsc(getcol(df,"PRF1"))), na.rm=TRUE)
    S <- ifn_score(df, excl=c(cands, CYT_SYMS))
    g <- getcol(df, cands); if(all(is.na(g))) next
    xs<-c(xs, abs(suppressWarnings(cor(g, S, method="spearman", use="complete.obs"))))
    ys<-c(ys, part(g, Y, S)) }
  data.frame(gene=sym, coreg=median(xs,na.rm=TRUE), partial=median(ys,na.rm=TRUE)) }
REAL <- dplyr::bind_rows(lapply(MARKER_SYMS, real_one))
write.csv(REAL, file.path(path_inter,"simulation_real_genes.csv"), row.names=FALSE)
cat("Done: 04_calibration_simulation.R\n")
