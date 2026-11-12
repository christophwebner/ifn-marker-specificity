# 05_proliferation_control.R
# Discrimination control in the four TCGA cohorts: the same test applied to
#   the interferon pathway (Reactome Interferon Signaling) with the thirteen markers and CYT, and
#   the E2F proliferation program (Hallmark E2F targets) with eight single-gene proliferation
#   markers and a G2/M mitotic read-out.
# The program score is the mean z-score of the program genes, without the tested gene if it is a
# member (leave-one-out) and without the read-out genes.
# Output (output/intermediate): proliferation_control.csv, proliferation_gene_sets.csv

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
source(file.path(root, "code", "gene_sets.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(janitor) })
set.seed(42)

# ---- Helpers (as in 02_specificity_test_tcga.R) ------------------------------
zsc  <- function(x){ s<-stats::sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s }
vv   <- function(df,ent){ c<-as.character(ent); if(!(c %in% colnames(df))) rep(NA_real_,nrow(df)) else df[[c]] }
sp   <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<20) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
part <- function(x,y,Z){ Z<-as.matrix(Z)
  ok<-is.finite(x)&is.finite(y)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  rx<-rank(x[ok]); ry<-rank(y[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  ex<-resid(lm(rx~RZ)); ey<-resid(lm(ry~RZ)); suppressWarnings(cor(ex,ey)) }
incR2 <- function(y,x,Z){ Z<-as.matrix(Z)
  ok<-is.finite(y)&is.finite(x)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  y<-rank(y[ok]); x<-rank(x[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  tryCatch(summary(lm(y~x+RZ))$r.squared - summary(lm(y~RZ))$r.squared, error=function(e) NA_real_) }
to_wide <- function(md){ md<-janitor::clean_names(as.data.frame(md))
  eg<-grep("entrez",colnames(md),value=TRUE)[1]; smp<-grep("^sample_id$|sampleid",colnames(md),value=TRUE)[1]; val<-grep("^value$",colnames(md),value=TRUE)[1]
  md %>% dplyr::transmute(sample_id=.data[[smp]], entrez=as.character(.data[[eg]]), value=.data[[val]]) %>%
    tidyr::pivot_wider(names_from=entrez, values_from=value) }

# ---- Programs, markers and read-outs (Entrez IDs) ---------------------------
.msigH <- tryCatch(as.data.frame(msigdbr::msigdbr(species="Homo sapiens", collection="H")),
                   error=function(e) as.data.frame(msigdbr::msigdbr(species="Homo sapiens", category="H")))
.ent <- grep("^ncbi_gene$|^entrez_gene$|ncbi_gene|entrez",colnames(.msigH),value=TRUE)[1]
.gsn <- grep("^gs_name$",colnames(.msigH),value=TRUE)[1]
hallmark <- function(setname){ v<-unique(suppressWarnings(as.integer(.msigH[[.ent]][.msigH[[.gsn]]==setname]))); v[!is.na(v)] }
PROG <- list(
  IFN    = ifn_entrez("primary"),
  PROLIF = hallmark("HALLMARK_E2F_TARGETS"))
stopifnot(all(vapply(PROG, length, 1L) > 50))
MARK <- list(
  IFN    = c(DDX58=23586, OAS1=4938, OASL=8638, GBP2=2634, IRF1=3659, USP18=11274, DDX60=55601,
             LGALS3BP=3959, ISG15=9636, IFI6=2537, IFI44=10561, IFI44L=10964, GBP1=2633),
  PROLIF = c(MKI67=4288, FOXM1=2305, MELK=9833, TK1=7083, CKS2=1164, BIRC5=332, ASPM=259266, EZH2=2146))
READ <- list(
  IFN    = c(GZMA=3001, PRF1=5551),
  PROLIF = c(CCNB1=891, CCNB2=9133, CDK1=983, CDC20=991, PLK1=5347, AURKA=6790, AURKB=9212,
             BUB1=699, BUB1B=701, TOP2A=7153, KIF11=3832, CENPA=1058, UBE2C=11065, TPX2=22974,
             NUSAP1=51203, KIF23=9493))
# read-out genes are removed from the program
for(pw in names(READ)){ READ[[pw]] <- READ[[pw]][!(READ[[pw]] %in% MARK[[pw]])]
                        PROG[[pw]] <- setdiff(PROG[[pw]], READ[[pw]]) }

.gs <- function(set, label, ents) data.frame(gene_set=set, description=label,
                                             gene_symbol=names(ents), entrez_id=unname(ents), stringsAsFactors=FALSE)
write.csv(rbind(.gs("G2M", "G2/M mitotic signature (read-out of the proliferation control)", READ$PROLIF),
                .gs("PROLIF_MARKER", "Single-gene proliferation markers (proliferation control)", MARK$PROLIF)),
          file.path(path_inter,"proliferation_gene_sets.csv"), row.names=FALSE)

# ---- TCGA PanCancer Atlas expression (cBioPortal, RSEM) --------------------
ALL <- unique(c(unlist(PROG), unlist(MARK), unlist(READ)))
STUDIES <- c(UCEC="ucec_tcga_pan_can_atlas_2018", OV="ov_tcga_pan_can_atlas_2018",
             CESC="cesc_tcga_pan_can_atlas_2018", BRCA="brca_tcga_pan_can_atlas_2018")
pull <- function(nm){ study<-STUDIES[nm]; prof<-paste0(study,"_rna_seq_v2_mrna")
  to_wide(cached_pull(paste0("tcga_",tolower(nm),"_proliferation_",genes_key(ALL)), function(){ cbio<-cBioPortalData::cBioPortal()
    cBioPortalData::getDataByGenes(cbio, studyId=study, genes=as.integer(unname(ALL)), by="entrezGeneId", molecularProfileIds=prof) })[[1]]) }
COH <- lapply(names(STUDIES), function(nm) pull(nm)); names(COH)<-names(STUDIES)

prog_score <- function(df, pw, exclude_ent=NA){
  ents<-PROG[[pw]]; keep <- if(is.na(exclude_ent)) rep(TRUE,length(ents)) else ents!=exclude_ent
  cols<-as.character(ents[keep]); cols<-cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_,nrow(df))); rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }
readout_vec <- function(df, pw){ cols<-as.character(READ[[pw]]); cols<-cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_,nrow(df))); rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }

# ---- Test per program, marker and cohort ------------------------------------
rows<-list(); k<-1
for(pw in names(PROG)){
  for(cn in names(COH)){ df<-COH[[cn]]; y<-readout_vec(df,pw)
    for(g in names(MARK[[pw]])){ ent<-MARK[[pw]][g]; x<-vv(df,ent); member<-ent %in% PROG[[pw]]
      sc<-prog_score(df, pw, exclude_ent=if(member) ent else NA)
      rows[[k]]<-data.frame(pathway=pw, cohort=cn, gene=g, type="published marker", member=member,
        rho_raw=round(sp(x,y),3), rho_partial=round(part(x,y,sc),3), unique_R2=round(incR2(y,x,sc),3)); k<-k+1 } }
}
R <- dplyr::bind_rows(rows)
write.csv(R, file.path(path_inter,"proliferation_control.csv"), row.names=FALSE)
cat("Done: 05_proliferation_control.R\n")
