# 02_specificity_test_tcga.R
# Specificity test in four TCGA cohorts: raw Spearman correlation and partial Spearman correlation
# after removal of the leave-one-out interferon pathway score (Reactome Interferon Signaling) for
# thirteen published markers and two controls against seven immune read-outs. Also computes the
# incremental R2, bootstrap confidence intervals (main read-outs), Meng's z-test with
# Benjamini-Hochberg correction, adjustment for leukocyte content and tumor purity, the
# classification of the markers and its sensitivity to the thresholds.
# Output (output/intermediate): readout_gene_sets.csv, tcga_statistics.csv,
#   tcga_classification.csv, threshold_sensitivity.csv

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
source(file.path(root, "code", "gene_sets.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(janitor) })
set.seed(42)
B_BOOT  <- 1500                 # bootstrap resamples
USE_ABS <- FALSE                # classification in the published (positive) direction

# ---- Genes (Entrez IDs) ----------------------------------------------------
CAND <- c(DDX58=23586, OAS1=4938, OASL=8638, GBP2=2634, IRF1=3659,
          USP18=11274, DDX60=55601, LGALS3BP=3959,
          ISG15=9636, IFI6=2537, IFI44=10561, IFI44L=10964, GBP1=2633)
MECH <- c(NLRC5=84166)
# Interferon pathway: Reactome (primary), Hallmark IFN-gamma and IFN-alpha (sensitivity)
REACTOME_IFN  <- ifn_entrez("primary")
HALLMARK_IFNG <- ifn_entrez("hallmark_ifng")
HALLMARK_IFNA <- ifn_entrez("hallmark_ifna")
IFN_DEFS <- list(
  reactome_ifn   = setNames(REACTOME_IFN,  REACTOME_IFN),
  hallmark_gamma = setNames(HALLMARK_IFNG, HALLMARK_IFNG),
  hallmark_alpha = setNames(HALLMARK_IFNA, HALLMARK_IFNA))
IFN_PRIMARY <- "reactome_ifn"
# Read-outs, positive control CD8A and leukocyte proxy PTPRC (CD45)
READ <- c(GZMA=3001, PRF1=5551, CD8A=925, CD274=29126,
          HLA_A=3105, HLA_B=3106, HLA_C=3107, B2M=567, PTPRC=5788,
          CD3D=915, CD3E=916, CD2=914, CD8B=926)
M2_GENES  <- c(CD163=9332, MRC1=4360, MSR1=4481, CD68=968, CSF1R=1436,
               VSIG4=11326, MS4A4A=51338, STAB1=23166, MARCO=8685, C1QC=714)
EXH_GENES <- c(PDCD1=5133, CTLA4=1493, LAG3=3902, HAVCR2=84868, TIGIT=201633,
               TOX=9760, ENTPD1=953, LAYN=143903)
APP_GENES <- c(TAP1=6890, TAP2=6891, TAPBP=6892, PSMB8=5696, PSMB9=5698, PSMB10=5699,
               ERAP1=51752, ERAP2=64167, CALR=811, CANX=821, PDIA3=2923,
               HLA_A=3105, HLA_B=3106, HLA_C=3107, B2M=567)
ALL  <- unique(c(CAND, MECH, unlist(IFN_DEFS), READ, M2_GENES, EXH_GENES, APP_GENES))

.gs <- function(set, label, ents) data.frame(gene_set=set, description=label,
                                             gene_symbol=sub("_","-",names(ents)), entrez_id=unname(ents), stringsAsFactors=FALSE)
GENESETS <- rbind(
  .gs("CYT",   "Cytolytic activity (main read-out)",            READ[c("GZMA","PRF1")]),
  .gs("EXH",   "T-cell exhaustion (main read-out)",             EXH_GENES),
  .gs("TCELL", "T-cell infiltration",                           READ[c("CD3D","CD3E","CD2","CD8B")]),
  .gs("MHC1",  "MHC class I",                                   READ[c("HLA_A","HLA_B","HLA_C","B2M")]),
  .gs("CD274", "PD-L1",                                         READ["CD274"]),
  .gs("M2",    "Myeloid / M2-like tumor-associated macrophages", M2_GENES),
  .gs("APP",   "Antigen processing and presentation",           APP_GENES),
  .gs("PTPRC", "Leukocyte proxy (confounder adjustment)",       READ["PTPRC"]),
  .gs("MARKER","Published single-gene markers",                 CAND),
  .gs("CTRL",  "Positive control (CD8A) and mechanism control (NLRC5)", c(READ["CD8A"], MECH["NLRC5"])))
write.csv(GENESETS, file.path(path_inter,"readout_gene_sets.csv"), row.names=FALSE)

# ---- TCGA PanCancer Atlas expression (cBioPortal, RSEM) --------------------
zsc <- function(x){ s<-stats::sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s }
to_wide <- function(md){ md<-janitor::clean_names(as.data.frame(md))
  eg<-grep("entrez",colnames(md),value=TRUE)[1]; smp<-grep("^sample_id$|sampleid",colnames(md),value=TRUE)[1]; val<-grep("^value$",colnames(md),value=TRUE)[1]
  md %>% dplyr::transmute(sample_id=.data[[smp]], entrez=as.character(.data[[eg]]), value=.data[[val]]) %>%
    tidyr::pivot_wider(names_from=entrez, values_from=value) }
STUDIES <- c(UCEC="ucec_tcga_pan_can_atlas_2018", OV="ov_tcga_pan_can_atlas_2018", CESC="cesc_tcga_pan_can_atlas_2018", BRCA="brca_tcga_pan_can_atlas_2018")
pull_targeted <- function(cache_name, study){ prof<-paste0(study,"_rna_seq_v2_mrna")
  cached_pull(cache_name, function(){ cbio<-cBioPortalData::cBioPortal()
    cBioPortalData::getDataByGenes(cbio, studyId=study, genes=as.integer(unname(ALL)), by="entrezGeneId", molecularProfileIds=prof) })[[1]] }
GENES_KEY <- genes_key(ALL)
uc <- to_wide(pull_targeted(paste0("tcga_ucec_main_",GENES_KEY),STUDIES["UCEC"]))
ov <- to_wide(pull_targeted(paste0("tcga_ov_main_",  GENES_KEY),STUDIES["OV"]))
ce <- to_wide(pull_targeted(paste0("tcga_cesc_main_",GENES_KEY),STUDIES["CESC"]))
br <- to_wide(pull_targeted(paste0("tcga_brca_main_",GENES_KEY),STUDIES["BRCA"]))
COH <- list(UCEC=uc, OV=ov, CESC=ce, BRCA=br)

# ---- Tumor purity (Aran et al. 2015, consensus CPE, ESTIMATE as fallback) --
ap <- suppressWarnings(readxl::read_excel(file.path(path_data,"aran_2015_purity.xlsx"), sheet="Supp Data 1", skip=2))
ap <- janitor::clean_names(as.data.frame(ap))
sidc <- grep("sample_id|^sample$", colnames(ap), value=TRUE)[1]
cpec <- grep("^cpe$", colnames(ap), value=TRUE); estc <- grep("^estimate$", colnames(ap), value=TRUE)
pv <- suppressWarnings(as.numeric(ap[[cpec[1]]]))
if(length(estc)){ ev<-suppressWarnings(as.numeric(ap[[estc[1]]])); pv<-ifelse(is.na(pv), ev, pv) }
PUR <- data.frame(key=substr(as.character(ap[[sidc]]),1,15), purity=pv, stringsAsFactors=FALSE)
PUR <- PUR[!is.na(PUR$key) & !is.na(PUR$purity) & !duplicated(PUR$key),]
attach_purity <- function(df){ df$purity <- PUR$purity[match(substr(df$sample_id,1,15), PUR$key)]; df }
COH <- lapply(COH, attach_purity)

need <- unique(c(CAND, MECH, READ))
for(nm in names(COH)){ miss <- names(need)[!(as.character(need) %in% colnames(COH[[nm]]))]
  if(length(miss)) warning(nm, ": genes missing from the download: ", paste(miss, collapse=", "), call.=FALSE) }

# ---- Statistics ------------------------------------------------------------
vv <- function(df,ent){ c<-as.character(ent); if(!(c %in% colnames(df))) rep(NA_real_,nrow(df)) else df[[c]] }
CYT  <- function(df) rowMeans(cbind(zsc(vv(df,READ["GZMA"])), zsc(vv(df,READ["PRF1"]))), na.rm=TRUE)
MHC1 <- function(df) rowMeans(cbind(zsc(vv(df,READ["HLA_A"])), zsc(vv(df,READ["HLA_B"])),
                                    zsc(vv(df,READ["HLA_C"])), zsc(vv(df,READ["B2M"]))), na.rm=TRUE)
TCELL <- function(df) rowMeans(cbind(zsc(vv(df,READ["CD3D"])), zsc(vv(df,READ["CD3E"])),
                                     zsc(vv(df,READ["CD2"])),  zsc(vv(df,READ["CD8B"]))), na.rm=TRUE)
# mean z-score of a gene list; exclude_ent removes the tested gene if it is a member (leave-one-out)
zmean_score <- function(df, ents, exclude_ent=NA){
  keep <- if(is.na(exclude_ent)) rep(TRUE,length(ents)) else ents!=exclude_ent
  cols <- as.character(ents[keep]); cols <- cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_, nrow(df)))
  rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }
M2  <- function(df, excl=NA) zmean_score(df, M2_GENES,  excl)
EXH <- function(df, excl=NA) zmean_score(df, EXH_GENES, excl)
APP <- function(df, excl=NA) zmean_score(df, APP_GENES, excl)
# Interferon pathway score without the genes in exclude_ent (tested gene and read-out genes)
ifn_score <- function(df, defn=IFN_PRIMARY, exclude_ent=NA){
  ents<-IFN_DEFS[[defn]]
  ex  <- suppressWarnings(as.integer(exclude_ent)); ex <- ex[!is.na(ex)]
  ents<-ents[!(ents %in% ex)]; cols<-as.character(ents); cols<-cols[cols %in% colnames(df)]
  if(!length(cols)) return(rep(NA_real_,nrow(df))); rowMeans(sapply(cols, function(c) zsc(df[[c]])), na.rm=TRUE) }
RO_GENES <- list(
  CYT   = READ[c("GZMA","PRF1")],
  EXH   = EXH_GENES,
  TCELL = READ[c("CD3D","CD3E","CD2","CD8B")],
  MHC1  = READ[c("HLA_A","HLA_B","HLA_C","B2M")],
  CD274 = READ["CD274"],
  M2    = M2_GENES,
  APP   = APP_GENES)
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<20) NA_real_ else suppressWarnings(cor(x[ok],y[ok],method="spearman")) }
# partial Spearman correlation (one or more covariates): correlation of rank residuals
part <- function(x,y,Z){ Z<-as.matrix(Z)
  ok<-is.finite(x)&is.finite(y)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  rx<-rank(x[ok]); ry<-rank(y[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  ex<-resid(lm(rx~RZ)); ey<-resid(lm(ry~RZ)); suppressWarnings(cor(ex,ey)) }
# incremental R2 of the gene beyond the covariates (ranks)
incR2 <- function(y,x,Z){ Z<-as.matrix(Z)
  ok<-is.finite(y)&is.finite(x)&apply(Z,1,function(r) all(is.finite(r))); if(sum(ok)<20) return(NA_real_)
  y<-rank(y[ok]); x<-rank(x[ok]); RZ<-apply(Z[ok,,drop=FALSE],2,rank)
  tryCatch(summary(lm(y~x+RZ))$r.squared - summary(lm(y~RZ))$r.squared, error=function(e) NA_real_) }
# Meng, Rosenthal and Rubin z-test for dependent overlapping correlations (marker vs. pathway score, same read-out)
meng_z <- function(r_xy, r_zy, r_xz, n){
  if(any(is.na(c(r_xy,r_zy,r_xz))) || n<10) return(c(z=NA_real_,p=NA_real_))
  rm2<-(r_xy^2+r_zy^2)/2; f<-(1-r_xz)/(2*(1-rm2)); f<-min(f,1); h<-(1-f*rm2)/(1-rm2)
  z<-(atanh(r_xy)-atanh(r_zy))*sqrt((n-3)/(2*(1-r_xz)*h)); c(z=z, p=2*pnorm(-abs(z))) }
# Bootstrap 95% CIs of the partial correlation, incremental R2 and reduction, computed on resampled
# ranks with the closed-form expressions
#   partial = (r_xy - r_xz r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))
#   incremental R2 = (r_xy^2 + r_yz^2 - 2 r_xy r_yz r_xz) / (1 - r_xz^2) - r_yz^2
boot_stats <- function(x,y,sc,B=B_BOOT){
  ok<-which(is.finite(x)&is.finite(y)&is.finite(sc)); n<-length(ok)
  if(n<30) return(c(rp_lo=NA,rp_hi=NA,uR2_lo=NA,uR2_hi=NA,drop_lo=NA,drop_hi=NA))
  rx<-rank(x[ok]); ry<-rank(y[ok]); rz<-rank(sc[ok])
  RP<-UR<-DR<-numeric(B)
  for(b in 1:B){ s<-sample.int(n,n,replace=TRUE)
    a<-rx[s]; c2<-ry[s]; d<-rz[s]
    rxy<-suppressWarnings(cor(a,c2)); rxz<-suppressWarnings(cor(a,d)); ryz<-suppressWarnings(cor(c2,d))
    den<-(1-rxz^2); denom<-sqrt(den*(1-ryz^2))
    rp<-if(is.finite(denom) && denom>1e-9) (rxy - rxz*ryz)/denom else NA_real_
    ur<-if(is.finite(den) && abs(den)>1e-9) max((rxy^2 + ryz^2 - 2*rxy*ryz*rxz)/den - ryz^2, 0) else NA_real_
    RP[b]<-rp; UR[b]<-ur; DR[b]<-rxy - rp }
  q<-function(v) as.numeric(quantile(v,c(.025,.975),na.rm=TRUE))
  c(rp_lo=q(RP)[1],rp_hi=q(RP)[2],uR2_lo=q(UR)[1],uR2_hi=q(UR)[2],drop_lo=q(DR)[1],drop_hi=q(DR)[2]) }

# ---- Specificity test --------------------------------------------------------
# excl_ro: tested gene, removed from the read-out score if it is a member
readout_vec <- function(df, r, excl_ro=NA){ if(r=="CYT") CYT(df) else if(r=="MHC1") MHC1(df) else if(r=="TCELL") TCELL(df) else if(r=="M2") M2(df,excl_ro) else if(r=="EXH") EXH(df,excl_ro) else if(r=="APP") APP(df,excl_ro) else vv(df, READ["CD274"]) }
run_one <- function(df, cn, r, gene, ent, typ, excl, do_boot){
  y<-readout_vec(df,r,excl_ro=ent); x<-vv(df,ent); sc<-ifn_score(df,IFN_PRIMARY,c(excl, RO_GENES[[r]])); cd45<-vv(df,READ["PTPRC"])
  pur<-if("purity" %in% colnames(df)) df$purity else rep(NA_real_,nrow(df))
  n<-sum(is.finite(x)&is.finite(y)&is.finite(sc))
  rr<-sp(x,y); rifn<-sp(sc,y); rxz<-sp(x,sc)
  rp<-part(x,y,sc); uR2<-incR2(y,x,sc); rp_leuko<-part(x,y,cbind(sc,cd45))
  rp_pur <- if(any(is.finite(pur))) part(x,y,cbind(sc,pur)) else NA_real_
  rp_pur_leuko <- if(any(is.finite(pur))) part(x,y,cbind(sc,pur,cd45)) else NA_real_
  mz<-meng_z(rr,rifn,rxz,n)
  ci<-if(do_boot) boot_stats(x,y,sc) else c(rp_lo=NA,rp_hi=NA,uR2_lo=NA,uR2_hi=NA,drop_lo=NA,drop_hi=NA)
  data.frame(cohort=cn, readout=r, gene=gene, type=typ, n=n,
             rho_raw=round(rr,3), rho_ifn=round(rifn,3),
             rho_partial=round(rp,3), rp_lo=round(ci["rp_lo"],3), rp_hi=round(ci["rp_hi"],3),
             unique_R2=round(uR2,3), uR2_lo=round(ci["uR2_lo"],3), uR2_hi=round(ci["uR2_hi"],3),
             drop=round(rr-rp,3), drop_lo=round(ci["drop_lo"],3), drop_hi=round(ci["drop_hi"],3),
             rho_partial_leuko=round(rp_leuko,3), rho_partial_purity=round(rp_pur,3),
             rho_partial_pur_leuko=round(rp_pur_leuko,3),
             meng_z=round(mz["z"],3), meng_p=signif(mz["p"],3),
             row.names=NULL) }

# Two main read-outs with bootstrap CIs, five further read-outs as point estimates
READOUTS      <- c("CYT","EXH","TCELL","MHC1","CD274","M2","APP")
MAIN_READOUTS <- c("CYT","EXH")
compute_R <- function(){
  res<-list(); k<-1
  for(cn in names(COH)){ df<-COH[[cn]]
    for(r in READOUTS){ db <- (r %in% MAIN_READOUTS)
      for(g in names(CAND)){ res[[k]]<-run_one(df,cn,r,g,CAND[g],"candidate (single-gene marker)",CAND[g],db); k<-k+1 } }
    # controls against both main read-outs, CD8A against T-cell infiltration, NLRC5 against MHC class I
    for(r in MAIN_READOUTS){
      res[[k]]<-run_one(df,cn,r,"CD8A", READ["CD8A"],  "positive control (CD8A)",   READ["CD8A"],  TRUE); k<-k+1
      res[[k]]<-run_one(df,cn,r,"NLRC5",MECH["NLRC5"],"mechanism control (NLRC5)",MECH["NLRC5"],TRUE); k<-k+1 }
    res[[k]]<-run_one(df,cn,"TCELL","CD8A", READ["CD8A"],"positive control (CD8A)",READ["CD8A"],  FALSE); k<-k+1
    res[[k]]<-run_one(df,cn,"MHC1", "NLRC5",MECH["NLRC5"],"mechanism control (NLRC5)",MECH["NLRC5"],TRUE); k<-k+1
  }
  R <- dplyr::bind_rows(res)
  R$meng_fdr <- NA_real_
  ci <- which(R$type=="candidate (single-gene marker)")
  R$meng_fdr[ci] <- signif(p.adjust(R$meng_p[ci], method="BH"),3)
  R$raw_real <- ifelse(abs(R$rho_raw)>=0.2, "yes", "weak(<0.2)")
  R }
R <- compute_R()

# ---- Classification ----------------------------------------------------------
# robust independent signal: median partial >= empirical threshold, partial >= 0.1 in every cohort,
#   median incremental R2 >= 0.05 and the same direction in at least 3/4 of the cohorts
# borderline: median partial >= empirical threshold with the same direction consistency
# proxy: all other markers
THRESHOLD <- as.numeric(read.csv(file.path(path_inter,"empirical_threshold.csv"))$threshold[1])
# partials: partial correlations of one gene across cohorts; use_abs = TRUE ignores the sign
tier_of <- function(partials, mur, use_abs=TRUE){
  p <- partials[is.finite(partials)]
  if(length(p) < 2 || !is.finite(mur)) return(NA_character_)
  meds <- median(p); sgn <- if(meds < 0) -1 else 1
  dir_frac <- mean(sign(p) == sgn | p == 0)
  med_mag  <- if(use_abs) median(abs(p)) else meds
  min_mag  <- if(use_abs) min(abs(p))    else min(p)
  bord_mag <- if(use_abs) median(abs(p)) else meds
  if(med_mag >= THRESHOLD && min_mag >= 0.1 && mur >= 0.05 && dir_frac >= 0.75) "robust independent signal"
  else if(bord_mag >= THRESHOLD && dir_frac >= 0.75)                    "inconsistent (borderline)"
  else                                                                  "IFN-pathway proxy" }
# TRUE if the class is unchanged whenever one cohort is left out
loo_tier_stable <- function(partials, urs, use_abs=TRUE){
  ok<-is.finite(partials); p<-partials[ok]; u<-urs[ok]; n<-length(p); if(n < 3) return(NA)
  base <- tier_of(p, median(u,na.rm=TRUE), use_abs)
  all(vapply(seq_len(n), function(i) identical(tier_of(p[-i], median(u[-i],na.rm=TRUE), use_abs), base), logical(1))) }
cls <- R %>% dplyr::filter(readout=="CYT", type=="candidate (single-gene marker)") %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(med_raw=round(median(rho_raw,na.rm=TRUE),3),
                   med_partial=round(median(rho_partial,na.rm=TRUE),3),
                   min_partial=round(min(rho_partial,na.rm=TRUE),3),
                   max_partial=round(max(rho_partial,na.rm=TRUE),3),
                   med_abs=round(median(abs(rho_partial),na.rm=TRUE),3),
                   dir_frac=round(mean(sign(rho_partial)==sign(median(rho_partial,na.rm=TRUE))),2),
                   med_uniqueR2=round(median(unique_R2,na.rm=TRUE),3),
                   med_partial_leuko=round(median(rho_partial_leuko,na.rm=TRUE),3),
                   med_partial_purity=round(median(rho_partial_purity,na.rm=TRUE),3),
                   direction    = ifelse(median(rho_partial,na.rm=TRUE)<0,"negative","positive"),
                   class_signed = tier_of(rho_partial, median(unique_R2,na.rm=TRUE), FALSE),
                   class_abs    = tier_of(rho_partial, median(unique_R2,na.rm=TRUE), TRUE),
                   loo_stable   = loo_tier_stable(rho_partial, unique_R2, USE_ABS),
                   .groups="drop") %>%
  dplyr::mutate(class = if(USE_ABS) class_abs else class_signed)
R <- R %>% dplyr::left_join(cls %>% dplyr::select(gene,class), by="gene") %>%
  dplyr::mutate(class=dplyr::case_when(
    type=="positive control (CD8A)"   ~ "positive control (CD8A)",
    type=="mechanism control (NLRC5)" ~ "mechanism control (NLRC5)", TRUE ~ class))
write.csv(R, file.path(path_inter,"tcga_statistics.csv"), row.names=FALSE)

# Classification of the thirteen markers on all seven read-outs (signed and absolute)
tier_tab <- R %>% dplyr::filter(type=="candidate (single-gene marker)") %>%
  dplyr::group_by(gene, readout) %>%
  dplyr::summarise(med_p =round(median(rho_partial,na.rm=TRUE),3),
                   med_ap=round(median(abs(rho_partial),na.rm=TRUE),3),
                   dir_frac=round(mean(sign(rho_partial)==sign(median(rho_partial,na.rm=TRUE))),2),
                   med_ur=round(median(unique_R2,na.rm=TRUE),3),
                   direction   = ifelse(median(rho_partial,na.rm=TRUE)<0,"-","+"),
                   tier_abs    = tier_of(rho_partial, median(unique_R2,na.rm=TRUE), TRUE),
                   tier_signed = tier_of(rho_partial, median(unique_R2,na.rm=TRUE), FALSE),
                   loo_stable  = loo_tier_stable(rho_partial, unique_R2, USE_ABS),
                   .groups="drop")
write.csv(tier_tab, file.path(path_inter,"tcga_classification.csv"), row.names=FALSE)

# ---- Threshold sensitivity (CYT) ---------------------------------------------
grid <- expand.grid(cut_partial=c(0.10,0.15,0.20,0.25,0.30), cut_uR2=c(0.02,0.05,0.10))
sens <- do.call(rbind, lapply(1:nrow(grid), function(i){
  retains <- cls %>% dplyr::filter(med_partial>=grid$cut_partial[i] & med_uniqueR2>=grid$cut_uR2[i]) %>% dplyr::pull(gene)
  data.frame(cut_partial=grid$cut_partial[i], cut_uR2=grid$cut_uR2[i],
             n_retains=length(retains), n_proxy=nrow(cls)-length(retains),
             retains=ifelse(length(retains)==0,"(none)",paste(retains,collapse=","))) }))
write.csv(sens, file.path(path_inter,"threshold_sensitivity.csv"), row.names=FALSE)
print(table(cls$class))
cat("Done: 02_specificity_test_tcga.R\n")
