# 06_figures_and_tables.R
# Builds all figures, tables, supplementary data and key numbers of the manuscript from the
# results of scripts 01 to 05 (output/intermediate).
# Output (output/paper, recreated on every run):
#   figures/   Fig1-Fig5, FigS1-FigS5 (PNG 300 dpi and PDF)
#   tables/    Table1, TableS2-TableS10 (CSV)
#   data/      DataS1-DataS3 (CSV)
#   key_numbers.csv, MANIFEST.txt

source(file.path(Sys.getenv("IFN_SPEC_ROOT", "."), "code", "setup.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

out_root <- file.path(path_output, "paper")
if(dir.exists(out_root)) unlink(out_root, recursive=TRUE)
dir_fig <- file.path(out_root, "figures"); dir_tab <- file.path(out_root, "tables"); dir_dat <- file.path(out_root, "data")
for(d in c(dir_fig, dir_tab, dir_dat)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
manifest <- character(0)
note <- function(file, desc) manifest <<- c(manifest, sprintf("%-50s %s", file, desc))

# ---- Helpers ---------------------------------------------------------------
rd <- function(f){ p <- file.path(path_inter, f)
  if(!file.exists(p)) stop("Missing input file: ", f, " (run scripts 01 to 05 first)", call.=FALSE)
  read.csv(p, stringsAsFactors=FALSE, check.names=FALSE) }
save_fig <- function(p, name, w, h, desc){
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width=w, height=h)
  note(paste0("figures/", name, " (.png, .pdf)"), desc) }
save_csv <- function(df, name, desc, dir=dir_tab, sub="tables"){
  write.csv(df, file.path(dir, paste0(name, ".csv")), row.names=FALSE)
  note(paste0(sub, "/", name, ".csv"), desc) }
relab <- function(x){ x <- as.character(x); x[x %in% c("RIG-I","RIGI")] <- "DDX58"; x }   # one label for DDX58
r3 <- function(x) round(x, 3)

THRESHOLD <- as.numeric(rd("empirical_threshold.csv")$threshold[1])

# ---- Labels and style ------------------------------------------------------
RO_LAB  <- c(CYT="Cytolytic activity (CYT)", EXH="T-cell exhaustion",
             TCELL="T-cell infiltration", MHC1="MHC class I", CD274="PD-L1 (CD274)",
             M2="Myeloid / M2-like macrophages", APP="Antigen processing and presentation")
MAIN_RO <- c("CYT","EXH")
SUPP_RO <- c("TCELL","MHC1","CD274","M2","APP")
# Classification in the published (positive) direction (tier_signed)
TIER_LAB <- c("robust independent signal"="Retains signal", "inconsistent (borderline)"="Borderline",
              "IFN-pathway proxy"="Pathway proxy")
# Fill and outline colors that remain distinguishable in greyscale
FILL <- c("Retains signal"="#B2182B", "Borderline"="#F4A582", "Pathway proxy"="white", "Control"="#5AAE61")
LINE <- c("Retains signal"="#B2182B", "Borderline"="#B35806", "Pathway proxy"="#2166AC", "Control"="#1B7837")
tier_scales <- function() list(scale_fill_manual(values=FILL), scale_colour_manual(values=LINE))
IFN_LAB <- "the interferon signaling pathway"
th  <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
         legend.title=element_blank(), strip.background=element_rect(fill="grey92", colour=NA),
         strip.text=element_text(face="bold"))
# theme for figures set at print width (6.9 inch)
th_print <- theme_bw(base_size=8) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
         legend.title=element_blank(), legend.text=element_text(size=7),
         strip.background=element_rect(fill="grey92", colour=NA), strip.text=element_text(face="bold", size=7),
         axis.text=element_text(size=6.5))
thr_lines <- function() list(geom_vline(xintercept=0, colour="grey60", linewidth=0.3),
                             geom_vline(xintercept=c(-THRESHOLD, THRESHOLD), linetype="dashed", colour="grey45", linewidth=0.4))

# ---- TCGA results (02) -------------------------------------------------------
R   <- rd("tcga_statistics.csv");     R$gene   <- relab(R$gene)
CLS <- rd("tcga_classification.csv"); CLS$gene <- relab(CLS$gene)
R   <- R %>% dplyr::filter(readout %in% c(MAIN_RO, SUPP_RO)) %>%
  dplyr::mutate(role=ifelse(grepl("control", type), "control", "marker"))

# median over the four TCGA cohorts per read-out and gene, with the classification
MED <- R %>% dplyr::group_by(readout, gene, role) %>%
  dplyr::summarise(n_cohorts=sum(is.finite(rho_partial)),
            raw=median(rho_raw, na.rm=TRUE), partial=median(rho_partial, na.rm=TRUE),
            partial_min=min(rho_partial, na.rm=TRUE), partial_max=max(rho_partial, na.rm=TRUE),
            unique_R2=median(unique_R2, na.rm=TRUE), .groups="drop") %>%
  dplyr::left_join(CLS %>% dplyr::select(gene, readout, tier_signed, loo_stable), by=c("gene","readout")) %>%
  dplyr::mutate(tier=ifelse(role=="control", "Control", unname(TIER_LAB[tier_signed])),
         direction=ifelse(partial < 0, "negative", "positive"))

# gene order: markers by median partial correlation with CYT, controls on top
ord_mark    <- MED %>% dplyr::filter(readout=="CYT", role=="marker") %>% dplyr::arrange(partial) %>% dplyr::pull(gene)
GENE_LEVELS <- c(ord_mark, "NLRC5", "CD8A")
fac <- function(d) d %>% dplyr::mutate(gene=factor(gene, levels=GENE_LEVELS))

# ---- Figure 1: specificity test on the two main read-outs -------------------
#   A raw vs. partial correlation (median over cohorts; grey dots: single cohorts)
#   B incremental R2 beyond the interferon pathway score
f1  <- MED %>% dplyr::filter(readout %in% MAIN_RO) %>% fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[MAIN_RO]))
f1c <- R   %>% dplyr::filter(readout %in% MAIN_RO) %>% fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[MAIN_RO]))
p1a <- ggplot(f1) + thr_lines() +
  geom_segment(aes(y=gene, yend=gene, x=raw, xend=partial), colour="grey75", linewidth=0.8) +
  geom_point(data=f1c, aes(x=rho_partial, y=gene), colour="grey55", size=0.9, alpha=0.7,
             position=position_jitter(height=0.12, width=0, seed=1)) +
  geom_point(aes(raw, gene), shape=16, colour="#22314E", size=2.2) +
  geom_point(aes(partial, gene, fill=tier, colour=tier), shape=24, size=2.8, stroke=0.8) +
  facet_wrap(~ro, nrow=1) + tier_scales() +
  labs(x=paste0("Spearman correlation with the read-out (circle: raw; triangle: partial, adjusted for ", IFN_LAB, ")"), y=NULL) + th + theme(axis.text.y=element_text(face="italic"))
p1b <- ggplot(f1, aes(unique_R2, gene, fill=tier, colour=tier)) + geom_col(width=0.65, linewidth=0.4) +
  geom_point(data=f1c, aes(x=unique_R2, y=gene), inherit.aes=FALSE, colour="grey30", size=0.8, alpha=0.7) +
  facet_wrap(~ro, nrow=1) + tier_scales() +
  labs(x=paste0("Incremental R-squared beyond ", IFN_LAB), y=NULL) + th + theme(legend.position="none", axis.text.y=element_text(face="italic"))
save_fig(p1a / p1b + plot_layout(heights=c(1.15, 1)) + plot_annotation(tag_levels="A", tag_prefix="(", tag_suffix=")"),
         "Fig1_specificity_test", 11, 11.5,
         "Specificity test, main read-outs CYT and T-cell exhaustion (A partial correlation, B incremental R2)")

# ---- Figure 2: each TCGA tumor type separately ------------------------------
COH_LAB <- c(UCEC="Endometrial (UCEC)", OV="Ovarian (OV)", CESC="Cervical (CESC)", BRCA="Breast (BRCA)")
n_coh   <- R %>% dplyr::filter(readout %in% MAIN_RO) %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(n=max(n, na.rm=TRUE), .groups="drop")
coh_lab <- setNames(sprintf("%s\nn = %d", COH_LAB[n_coh$cohort], n_coh$n), n_coh$cohort)
fT <- R %>% dplyr::filter(readout %in% MAIN_RO) %>%
  dplyr::left_join(MED %>% dplyr::filter(readout %in% MAIN_RO) %>% dplyr::select(readout, gene, tier), by=c("readout","gene")) %>%
  fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[MAIN_RO]),
                          cohort_lab=factor(coh_lab[cohort], levels=unname(coh_lab[names(COH_LAB)])))
figT <- ggplot(fT) + thr_lines() +
  geom_segment(aes(y=gene, yend=gene, x=rho_raw, xend=rho_partial), colour="grey75", linewidth=0.5) +
  geom_errorbar(aes(y=gene, xmin=rp_lo, xmax=rp_hi, colour=tier), width=0.35, linewidth=0.35, orientation="y",
                na.rm=TRUE, show.legend=FALSE) +
  geom_point(aes(rho_raw, gene), shape=16, colour="#22314E", size=1.2) +
  geom_point(aes(rho_partial, gene, fill=tier, colour=tier), shape=24, size=1.8, stroke=0.6) +
  facet_grid(ro ~ cohort_lab, labeller=labeller(ro=label_wrap_gen(14))) + tier_scales() +
  scale_x_continuous(breaks=c(-0.4, 0, 0.4, 0.8)) +
  labs(x="Spearman correlation with the read-out (circle: raw; triangle: partial with 95% bootstrap CI)", y=NULL) + th_print + theme(axis.text.y=element_text(face="italic"))
save_fig(figT, "Fig2_tumor_types", 6.9, 7.4,
         "Specificity test in each TCGA tumor type separately (endometrial, ovarian, cervical, breast); read-outs CYT and T-cell exhaustion")

# ---- Figure 3: calibration simulation with the real genes (04) --------------
SIM  <- rd("simulation_grid.csv")
AGG  <- SIM %>% dplyr::group_by(a, b) %>%
  dplyr::summarise(median_partial=median(median_partial, na.rm=TRUE), coreg=median(coreg, na.rm=TRUE),
            pass_rate=mean(pass_rate, na.rm=TRUE), .groups="drop") %>% dplyr::mutate(b_lab=paste0("b = ", b))
REAL <- rd("simulation_real_genes.csv") %>% dplyr::mutate(gene=relab(gene)) %>% dplyr::select(gene, coreg, partial) %>%
  dplyr::left_join(MED %>% dplyr::filter(readout=="CYT") %>% dplyr::select(gene, tier), by="gene")
lab_pos <- AGG %>% dplyr::group_by(b_lab) %>% dplyr::slice_min(coreg, n=1, with_ties=FALSE) %>% dplyr::ungroup()
fig2 <- ggplot() +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.3) +
  geom_hline(yintercept=c(-THRESHOLD, THRESHOLD), linetype="dashed", colour="grey45", linewidth=0.4) +
  geom_line(data=AGG, aes(coreg, median_partial, group=b_lab), colour="grey65", linewidth=0.8) +
  geom_point(data=AGG, aes(coreg, median_partial), colour="grey65", size=1.4) +
  geom_text(data=lab_pos, aes(coreg, median_partial, label=b_lab), hjust=1.15, size=3, colour="grey35") +
  geom_point(data=REAL, aes(coreg, partial, fill=tier, colour=tier), shape=21, size=2.8, stroke=0.8) +
  ggrepel::geom_text_repel(data=REAL, aes(coreg, partial, label=gene, colour=tier), size=3, fontface="italic", seed=1,
                           max.overlaps=30, show.legend=FALSE) +
  tier_scales() + scale_x_continuous(expand=expansion(mult=c(0.12, 0.05))) +
  labs(x=paste0("Co-regulation with ", IFN_LAB, " (|Spearman correlation| with the pathway score)"),
       y="Partial correlation with cytolytic activity") + th
save_fig(fig2, "Fig3_calibration_simulation", 10, 6.8,
         "Semi-synthetic calibration: synthetic genes with known independent signal b (grey) and real genes")

# ---- Figure 4: the whole interferon signaling pathway (01) -----------------
# Pathway genes are counted only if their raw correlation reaches the threshold (source data: Table S10).
PW_RAW_MIN <- THRESHOLD
BGW  <- rd("background_correlations.csv") %>% dplyr::mutate(gene=relab(gene))
n_pw <- sum(BGW$role=="IFN-programme gene"); n_mt <- sum(BGW$role=="matched non-IFN gene")
GRP_MT  <- sprintf("Expression-matched non-interferon genes (n = %d)", n_mt)
GRP_IFN <- sprintf("Interferon signaling pathway genes (Reactome, n = %d)", n_pw)
f4cloud <- BGW %>% dplyr::filter(role %in% c("IFN-programme gene", "matched non-IFN gene")) %>%
  dplyr::mutate(grp=factor(ifelse(role=="IFN-programme gene", GRP_IFN, GRP_MT), levels=c(GRP_MT, GRP_IFN))) %>%
  dplyr::arrange(grp)
f4mk <- BGW %>% dplyr::filter(role=="published marker") %>%
  dplyr::left_join(MED %>% dplyr::filter(readout=="CYT", role=="marker") %>% dplyr::select(gene, tier), by="gene")
lim4 <- range(c(f4cloud$raw, f4cloud$partial), na.rm=TRUE) + c(-0.02, 0.02)
th4 <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.title=element_blank(),
         legend.text=element_text(size=8), legend.key.size=unit(10, "pt"), legend.key.spacing.y=unit(5, "pt"))
GRP_MT_W  <- sub(" \\(", "\n(", GRP_MT); GRP_IFN_W <- sub(" \\(", "\n(", GRP_IFN)
p4a <- ggplot(f4cloud, aes(raw, partial)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.3) + geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
  geom_hline(yintercept=c(-THRESHOLD, THRESHOLD), linetype="dashed", colour="grey45", linewidth=0.35) +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey20", linewidth=0.35) +
  geom_point(aes(colour=grp), size=1, alpha=0.75, stroke=0) +
  geom_point(data=f4mk, aes(raw, partial, fill=tier), shape=21, colour="black", size=2.4, stroke=0.5) +
  ggrepel::geom_text_repel(data=f4mk, aes(raw, partial, label=gene), size=2.6, fontface="bold.italic", colour="grey10",
                           seed=1, max.overlaps=Inf, min.segment.length=0, box.padding=0.3, point.padding=0.15,
                           segment.colour="grey40") +
  scale_colour_manual(values=setNames(c("grey78", "#2F5F98"), c(GRP_MT, GRP_IFN)), breaks=c(GRP_MT, GRP_IFN),
                      labels=c(GRP_MT_W, GRP_IFN_W)) +
  scale_fill_manual(values=FILL, labels=function(x) paste0("Published marker,\n", tolower(x))) +
  coord_equal(xlim=lim4, ylim=lim4) +
  labs(x="Raw correlation with cytolytic activity", y="Partial correlation after adjustment\nfor the pathway score") + th4 +
  theme(legend.position="right", legend.box="vertical", legend.spacing.y=unit(8, "pt")) +
  guides(colour=guide_legend(ncol=1, order=1, override.aes=list(size=2.2, alpha=1)), fill=guide_legend(ncol=1, order=2))
save_fig(p4a, "Fig4_pathway_membership", 6.3, 4.6,
         "All Reactome interferon genes vs. expression-matched non-interferon genes, raw vs. partial correlation with CYT; published markers highlighted")
f4b <- BGW %>% dplyr::filter(role=="IFN-programme gene", raw >= PW_RAW_MIN) %>% dplyr::arrange(dplyr::desc(partial))
n_below <- sum(f4b$partial < THRESHOLD); n_above <- sum(f4b$partial >= THRESHOLD)

# ---- Figure 5: external validation (03) -------------------------------------
V  <- rd("external_validation.csv"); V$gene <- relab(V$gene)
V3 <- V %>% dplyr::filter(readout %in% MAIN_RO) %>%
  dplyr::left_join(MED %>% dplyr::filter(readout %in% MAIN_RO) %>% dplyr::select(readout, gene, tier), by=c("readout","gene")) %>%
  fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[MAIN_RO]), cohort_lab=factor(label, levels=unique(label)))
fig3 <- ggplot(V3) + thr_lines() +
  geom_segment(aes(y=gene, yend=gene, x=rho_raw, xend=rho_partial), colour="grey75", linewidth=0.5) +
  geom_errorbar(aes(y=gene, xmin=rp_lo, xmax=rp_hi, colour=tier), width=0.35, linewidth=0.35, orientation="y",
                na.rm=TRUE, show.legend=FALSE) +
  geom_point(aes(rho_raw, gene), shape=16, colour="#22314E", size=1.2) +
  geom_point(aes(rho_partial, gene, fill=tier, colour=tier), shape=24, size=1.8, stroke=0.6) +
  facet_grid(ro ~ cohort_lab, labeller=labeller(ro=label_wrap_gen(14), cohort_lab=label_wrap_gen(20))) + tier_scales() +
  scale_x_continuous(breaks=c(-0.4, 0, 0.4, 0.8)) +
  labs(x="Spearman correlation with the read-out (circle: raw; triangle: partial with 95% bootstrap CI)", y=NULL) + th_print + theme(axis.text.y=element_text(face="italic"))
save_fig(fig3, "Fig5_external_validation", 6.9, 7.4,
         "External validation: METABRIC, CPTAC-UCEC, GSE9891; read-outs CYT and T-cell exhaustion")

# ---- Figure S1: additional read-outs ----------------------------------------
s1  <- MED %>% dplyr::filter(readout %in% SUPP_RO) %>% fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[SUPP_RO]))
s1c <- R   %>% dplyr::filter(readout %in% SUPP_RO) %>% fac() %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[SUPP_RO]))
figS1 <- ggplot(s1) + thr_lines() +
  geom_segment(aes(y=gene, yend=gene, x=raw, xend=partial), colour="grey75", linewidth=0.8) +
  geom_point(data=s1c, aes(x=rho_partial, y=gene), colour="grey55", size=0.8, alpha=0.7,
             position=position_jitter(height=0.12, width=0, seed=1)) +
  geom_point(aes(raw, gene), shape=16, colour="#22314E", size=2) +
  geom_point(aes(partial, gene, fill=tier, colour=tier), shape=24, size=2.5, stroke=0.7) +
  facet_wrap(~ro, nrow=2) + tier_scales() +
  labs(x="Spearman correlation with the read-out (circle: raw; triangle: partial)", y=NULL) + th + theme(axis.text.y=element_text(face="italic"))
save_fig(figS1, "FigS1_additional_readouts", 12, 9,
         "Specificity test on five additional read-outs")

# ---- Figure S2: constructions of the interferon pathway score (01) ----------
SC <- rd("score_construction.csv"); SC$gene <- relab(SC$gene)
SC_RO <- c("CYT","MHC1","APP")
SC <- SC %>% dplyr::filter(readout %in% SC_RO) %>% dplyr::mutate(ro=factor(RO_LAB[readout], levels=RO_LAB[SC_RO]))
SC_LAB <- c(primary_meanz      ="Reactome Interferon Signaling, mean z-score (primary)",
            primary_pc1        ="Reactome Interferon Signaling, first principal component",
            primary_gsva       ="Reactome Interferon Signaling, GSVA",
            primary_ssgsea     ="Reactome Interferon Signaling, ssGSEA",
            hallmark_ifng_meanz="Hallmark IFN-gamma response, mean z-score",
            hallmark_ifna_meanz="Hallmark IFN-alpha response, mean z-score")
s2 <- SC %>% tidyr::pivot_longer(dplyr::all_of(names(SC_LAB)), names_to="score", values_to="partial") %>%
  dplyr::mutate(score=factor(SC_LAB[score], levels=SC_LAB)) %>% fac()
s2raw <- SC %>% fac()
figS2 <- ggplot(s2, aes(partial, gene, colour=score, shape=score)) + thr_lines() +
  geom_point(data=s2raw, aes(raw, gene), inherit.aes=FALSE, shape=1, size=2.6, colour="grey20") +
  geom_point(size=2.3, position=position_dodge(width=0.7)) +
  facet_wrap(~ro, nrow=1, scales="free_y") +
  scale_x_continuous(breaks=round(seq(-0.6, 1, 0.2), 1), labels=function(x) sprintf("%.1f", x)) +
  scale_colour_brewer(palette="Dark2") + scale_shape_manual(values=c(17,16,15,18,8,4)) +
  labs(x="Correlation with the read-out (median across cohorts; open circle: unadjusted; symbols: partial)", y=NULL) + th + theme(axis.text.y=element_text(face="italic")) +
  guides(colour=guide_legend(nrow=3), shape=guide_legend(nrow=3))
save_fig(figS2, "FigS2_score_construction", 16, 9,
         "Robustness to the interferon pathway score (Reactome: mean z, PC1, GSVA, ssGSEA; Hallmark IFN-gamma and IFN-alpha) on CYT, MHC class I and antigen presentation; controls CD8A and NLRC5")

# ---- Figure S3: adjustment for leukocyte content and tumor purity (02) ------
ADJ_LAB <- c(rho_partial="Interferon signaling pathway",
             rho_partial_leuko="+ leukocyte fraction (PTPRC)",
             rho_partial_purity="+ tumor purity",
             rho_partial_pur_leuko="+ tumor purity + leukocyte fraction")
s3 <- R %>% dplyr::filter(readout %in% MAIN_RO) %>% dplyr::select(readout, gene, cohort, dplyr::all_of(names(ADJ_LAB))) %>%
  tidyr::pivot_longer(dplyr::all_of(names(ADJ_LAB)), names_to="adj", values_to="partial") %>%
  dplyr::group_by(readout, gene, adj) %>% dplyr::summarise(partial=median(partial, na.rm=TRUE), .groups="drop") %>%
  dplyr::mutate(adj=factor(ADJ_LAB[adj], levels=ADJ_LAB), ro=factor(RO_LAB[readout], levels=RO_LAB[MAIN_RO])) %>% fac()
figS3 <- ggplot(s3, aes(partial, gene, colour=adj, shape=adj)) + thr_lines() +
  geom_point(size=2.2, position=position_dodge(width=0.6), na.rm=TRUE) +
  facet_wrap(~ro, nrow=1) + scale_colour_brewer(palette="Set1") + scale_shape_manual(values=c(17,16,15,18)) +
  labs(x="Partial correlation with the read-out (median across cohorts)", y=NULL) + th + theme(axis.text.y=element_text(face="italic")) +
  guides(colour=guide_legend(nrow=2), shape=guide_legend(nrow=2))
save_fig(figS3, "FigS3_confounder_adjustment", 11, 7.5,
         "Additional adjustment for leukocyte fraction and tumor purity")

# ---- Figure S4: derivation of the empirical threshold (01) ------------------
BG <- rd("background_correlations.csv"); BG$gene <- relab(BG$gene)
ROLE_LAB <- c("matched non-IFN gene"="Expression-matched non-interferon genes",
              "IFN-programme gene"="Interferon signaling pathway genes (Reactome)",
              "published marker"="Published single-gene markers")
bgd <- BG %>% dplyr::filter(role %in% names(ROLE_LAB)) %>% dplyr::mutate(role=factor(ROLE_LAB[role], levels=ROLE_LAB))
figS4 <- ggplot(bgd %>% dplyr::filter(role != ROLE_LAB[["published marker"]]), aes(partial, fill=role)) +
  geom_density(alpha=0.45, colour=NA) +
  geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
  geom_vline(xintercept=c(-THRESHOLD, THRESHOLD), linetype="dashed", colour="grey25", linewidth=0.5) +
  geom_rug(data=bgd %>% dplyr::filter(role == ROLE_LAB[["published marker"]]), aes(x=partial), inherit.aes=FALSE,
           colour="black", linewidth=0.8, length=unit(0.06, "npc")) +
  scale_fill_manual(values=c("grey78", "#2F5F98")) +
  labs(x=sprintf("Partial correlation with cytolytic activity (dashed: empirical threshold |partial| = %.2f; ticks: published markers)", THRESHOLD),
       y="Density") + th
save_fig(figS4, "FigS4_threshold_derivation", 10, 5.5,
         "Null distribution of expression-matched non-interferon genes and the empirical threshold")

# ---- Figure S5: interferon pathway vs. proliferation program (05) -----------
PG <- rd("proliferation_control.csv"); PG$gene <- relab(PG$gene)
PW_LAB <- c(IFN="Interferon signaling pathway, read-out cytolytic activity",
            PROLIF="E2F proliferation program, read-out G2/M mitotic signature")
s5 <- PG %>% dplyr::filter(type=="published marker", pathway %in% names(PW_LAB)) %>% dplyr::group_by(pathway, gene) %>%
  dplyr::summarise(raw=median(rho_raw, na.rm=TRUE), partial=median(rho_partial, na.rm=TRUE),
            unique_R2=median(unique_R2, na.rm=TRUE), .groups="drop") %>%
  dplyr::mutate(status=ifelse(partial < THRESHOLD, "Loses signal", "Retains signal"),
         pw=factor(PW_LAB[pathway], levels=PW_LAB), key=paste(pathway, gene, sep="__"))
s5$key <- factor(s5$key, levels=s5$key[order(s5$partial)])
figS5 <- ggplot(s5) + thr_lines() +
  geom_segment(aes(y=key, yend=key, x=raw, xend=partial), colour="grey75", linewidth=0.8) +
  geom_point(aes(raw, key), shape=16, colour="#22314E", size=2.2) +
  geom_point(aes(partial, key, fill=status, colour=status), shape=24, size=2.8, stroke=0.8) +
  facet_wrap(~pw, scales="free_y") + scale_y_discrete(labels=function(x) sub("^.*__", "", x)) +
  scale_fill_manual(values=c("Loses signal"="white", "Retains signal"=unname(FILL["Retains signal"]))) +
  scale_colour_manual(values=c("Loses signal"=unname(LINE["Pathway proxy"]), "Retains signal"=unname(LINE["Retains signal"]))) +
  labs(x="Spearman correlation with the read-out (circle: raw; triangle: partial, adjusted for the pathway or program)", y=NULL) + th + theme(axis.text.y=element_text(face="italic"))
save_fig(figS5, "FigS5_program_discrimination", 11, 5.5,
         "Discrimination control: interferon signaling pathway vs. proliferation program")

# ---- Tables ----------------------------------------------------------------
# Table 1: classification on the two main read-outs (median over the four TCGA cohorts)
t1cols <- c("raw","partial","partial_min","partial_max","unique_R2","direction","classification")
t1 <- MED %>% dplyr::filter(readout %in% MAIN_RO) %>%
  dplyr::transmute(gene, role, readout, raw=r3(raw), partial=r3(partial), partial_min=r3(partial_min),
            partial_max=r3(partial_max), unique_R2=r3(unique_R2), direction, classification=tier) %>%
  tidyr::pivot_wider(names_from=readout, values_from=dplyr::all_of(t1cols), names_glue="{readout}_{.value}")
t1 <- t1[order(match(t1$gene, rev(GENE_LEVELS))), c("gene","role", as.vector(t(outer(MAIN_RO, t1cols, paste, sep="_"))))]
save_csv(t1, "Table1_classification", "Classification of 13 markers and 2 controls on CYT and T-cell exhaustion")

# Table S2: gene sets of the read-outs, leukocyte proxy, markers, controls and proliferation control
save_csv(rbind(rd("readout_gene_sets.csv"), rd("proliferation_gene_sets.csv")), "TableS2_gene_sets",
         "Gene sets of all read-outs, the leukocyte proxy, the markers, the controls and the proliferation control")

# Table S3: per-cohort statistics (TCGA)
t3 <- R %>% dplyr::transmute(cohort, readout, readout_label=RO_LAB[readout], gene, role, n,
                      rho_raw, rho_ifn_pathway=rho_ifn, rho_partial, rho_partial_ci_lo=rp_lo, rho_partial_ci_hi=rp_hi,
                      unique_R2, unique_R2_ci_lo=uR2_lo, unique_R2_ci_hi=uR2_hi,
                      rho_partial_adj_leukocytes=rho_partial_leuko, rho_partial_adj_purity=rho_partial_purity,
                      rho_partial_adj_purity_leukocytes=rho_partial_pur_leuko, meng_z, meng_p, meng_fdr) %>%
  dplyr::arrange(match(readout, c(MAIN_RO, SUPP_RO)), cohort, gene)
save_csv(t3, "TableS3_full_statistics", "Per-cohort statistics for all genes and read-outs (TCGA)")

# Table S4: classification on all seven read-outs, in the published direction and ignoring the sign
t4 <- CLS %>% dplyr::filter(readout %in% c(MAIN_RO, SUPP_RO)) %>%
  dplyr::transmute(gene, readout, readout_label=RO_LAB[readout], median_partial=med_p, median_abs_partial=med_ap,
            direction_consistency=dir_frac, median_unique_R2=med_ur, direction,
            classification=unname(TIER_LAB[tier_signed]), classification_ignoring_sign=unname(TIER_LAB[tier_abs]),
            stable_leave_one_cohort_out=loo_stable) %>%
  dplyr::arrange(match(readout, c(MAIN_RO, SUPP_RO)), match(gene, rev(GENE_LEVELS)))
save_csv(t4, "TableS4_classification_all_readouts", "Classification of the 13 markers on all seven read-outs")

# Table S5: sensitivity of the classification to the thresholds
save_csv(rd("threshold_sensitivity.csv"), "TableS5_threshold_sensitivity", "Number of markers retaining signal across threshold combinations (CYT)")

# Table S6: constructions of the interferon pathway score
t6 <- SC %>% dplyr::transmute(readout, readout_label=RO_LAB[readout], gene, role, raw, reactome_mean_z=primary_meanz,
                    reactome_pc1=primary_pc1, reactome_gsva=primary_gsva,
                    reactome_ssgsea=primary_ssgsea, hallmark_ifng_mean_z=hallmark_ifng_meanz,
                    hallmark_ifna_mean_z=hallmark_ifna_meanz) %>%
  dplyr::arrange(match(readout, SC_RO), match(gene, rev(GENE_LEVELS)))
save_csv(t6, "TableS6_score_construction", "Median raw and partial correlation under six constructions of the IFN pathway score (CYT, MHC class I, antigen presentation; controls CD8A and NLRC5)")

# Table S7: external validation
t7 <- V %>% dplyr::filter(readout %in% MAIN_RO | (gene=="NLRC5" & readout=="MHC1")) %>%
  dplyr::transmute(cohort=label, readout, readout_label=RO_LAB[readout], gene, type, n, rho_raw, rho_partial,
            rho_partial_ci_lo=rp_lo, rho_partial_ci_hi=rp_hi, unique_R2)
save_csv(t7, "TableS7_external_validation", "External validation statistics (METABRIC, CPTAC-UCEC, GSE9891)")

# Table S8: calibration simulation
t8 <- AGG %>% dplyr::transmute(coregulation_loading_a=a, independent_signal_loading_b=b, measured_coregulation=r3(coreg),
                        median_partial=r3(median_partial), pass_rate=r3(pass_rate))
save_csv(t8, "TableS8_simulation_calibration", "Simulation grid: median partial correlation and pass rate per (a, b)")

# Table S9: interferon pathway vs. proliferation program
t9 <- s5 %>% dplyr::transmute(program=PW_LAB[pathway], gene, raw=r3(raw), partial=r3(partial), unique_R2=r3(unique_R2), status) %>%
  dplyr::arrange(program, partial)
save_csv(t9, "TableS9_program_discrimination", "Interferon signaling pathway vs. proliferation program: median raw and partial correlations")

# Table S10: every Reactome interferon gene and its expression-matched partner (source data of Figure 4)
t10 <- BGW %>% dplyr::filter(role %in% c("IFN-programme gene", "matched non-IFN gene")) %>%
  dplyr::transmute(gene, group=ifelse(role=="IFN-programme gene", "Interferon signaling pathway gene (Reactome)",
                                      "Expression-matched non-interferon gene"),
                   published_marker=gene %in% f4mk$gene, raw=r3(raw), partial=r3(partial), reduction=r3(raw - partial),
                   raw_at_least_min=raw >= PW_RAW_MIN,
                   partial_vs_threshold=ifelse(abs(partial) >= THRESHOLD, "above threshold", "below threshold")) %>%
  dplyr::arrange(group, dplyr::desc(partial))
save_csv(t10, "TableS10_pathway_genes", sprintf("Raw and partial correlation with CYT for all Reactome interferon genes and their expression-matched partners (source data of Figure 4; threshold %.2f, raw >= %.2f flagged)", THRESHOLD, PW_RAW_MIN))

# ---- Supplementary data ------------------------------------------------------
save_csv(rd("ifn_gene_sets.csv"), "DataS1_ifn_gene_sets",
         "Interferon gene sets: Reactome Interferon Signaling (primary), Hallmark IFN-gamma and IFN-alpha response", dir_dat, "data")
save_csv(read.csv(file.path(path_data, "matched_gene_pairs.csv"), check.names=FALSE), "DataS2_matched_gene_pairs",
         "Interferon signaling pathway genes and their expression-matched non-interferon partners", dir_dat, "data")
save_csv(read.csv(file.path(path_data, "background_gene_pool.csv"), check.names=FALSE), "DataS3_background_gene_pool",
         "Fixed pool of non-interferon background genes used for matching", dir_dat, "data")

# ---- Key numbers for the manuscript text -------------------------------------
tier_n <- function(r) { x <- MED %>% dplyr::filter(readout==r, role=="marker"); vapply(unname(TIER_LAB), function(t) sum(x$tier==t, na.rm=TRUE), numeric(1)) }
red_pw <- BGW$reduction[BGW$role=="IFN-programme gene"]; red_mt <- BGW$reduction[BGW$role=="matched non-IFN gene"]
kn <- rbind(
  data.frame(item="empirical_threshold_abs_partial", value=r3(THRESHOLD)),
  as.data.frame(R %>% dplyr::filter(readout=="CYT", role=="marker") %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(value=max(n, na.rm=TRUE), .groups="drop") %>% dplyr::transmute(item=paste0("n_tumors_", cohort), value)),
  data.frame(item=paste0("CYT_n_", gsub(" ", "_", unname(TIER_LAB))), value=unname(tier_n("CYT"))),
  data.frame(item=paste0("EXH_n_", gsub(" ", "_", unname(TIER_LAB))), value=unname(tier_n("EXH"))),
  data.frame(item=c("pathway_genes_total", "pathway_genes_raw_ge_min", "pathway_genes_below_threshold",
                    "pathway_genes_above_threshold", "matched_genes_raw_ge_min", "matched_genes_raw_ge_min_below_threshold"),
             value=c(n_pw, nrow(f4b), n_below, n_above,
                     sum(BGW$role=="matched non-IFN gene" & BGW$raw >= PW_RAW_MIN),
                     sum(BGW$role=="matched non-IFN gene" & BGW$raw >= PW_RAW_MIN & BGW$partial < THRESHOLD))),
  data.frame(item=c("median_reduction_pathway_genes_raw_ge_min", "median_reduction_matched_genes_raw_ge_min",
                    "wilcoxon_p_reduction_pathway_vs_matched_all_genes"),
             value=c(r3(median(f4b$reduction, na.rm=TRUE)),
                     r3(median(BGW$reduction[BGW$role=="matched non-IFN gene" & BGW$raw >= PW_RAW_MIN], na.rm=TRUE)),
                     signif(wilcox.test(red_pw, red_mt, exact=FALSE)$p.value, 3))),
  data.frame(item="simulation_false_positive_rate_b0", value=r3(mean(AGG$pass_rate[AGG$b==0]))),
  data.frame(item="simulation_max_span_over_a", value=r3(max(tapply(AGG$median_partial, AGG$b, function(v) max(v)-min(v))))),
  as.data.frame(V %>% dplyr::group_by(label) %>% dplyr::summarise(value=max(n, na.rm=TRUE), .groups="drop") %>%
    dplyr::transmute(item=paste0("n_external_", gsub("[^A-Za-z0-9]+", "_", label)), value)),
  as.data.frame(s5 %>% dplyr::group_by(pathway) %>% dplyr::summarise(value=r3(mean(partial < THRESHOLD)), .groups="drop") %>%
    dplyr::transmute(item=paste0("fraction_losing_signal_", pathway), value)))
write.csv(kn, file.path(out_root, "key_numbers.csv"), row.names=FALSE)
note("key_numbers.csv", "Reference values for the manuscript text (not a figure or table)")

writeLines(c(paste("Paper outputs generated", format(Sys.time(), "%Y-%m-%d %H:%M")), "", manifest),
           file.path(out_root, "MANIFEST.txt"))
print(kn, row.names=FALSE)
cat("Done: 06_figures_and_tables.R, outputs in", out_root, "\n")
