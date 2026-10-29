# =============================================================================
# analysis.R — differential expression (sleuth), GSEA (fgsea + MSigDB),
# publication figures and tables for the GSE79013 reanalysis.
# Called by run_all.sh; can also be run alone from this folder:  Rscript analysis.R
# =============================================================================
suppressPackageStartupMessages({
  library(sleuth); library(data.table); library(fgsea); library(msigdbr)
  library(ggplot2); library(ggrepel); library(patchwork); library(pheatmap)
  library(jsonlite)
})
set.seed(42)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

FDR_GENE  <- 0.05      # gene-level q-value cutoff
FDR_PATH  <- 0.05      # pathway padj cutoff used for figures (all results kept in Table S2)
KEY_GENES <- c("Nppb", "Acta1", "Abra", "Ech1", "Wee1", "Cacng6")   # Figure 3B (if present)
# Primary test: sleuth Wald test on the treatment coefficient (as in the original analysis);
# the likelihood-ratio test (LRT) is reported alongside as a more conservative sensitivity analysis.
# Pathways shown in Figure 2B/C (those discussed in the manuscript); falls back to top-3 if not significant.
SHOW_UP   <- c("REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES", "REACTOME_COLLAGEN_FORMATION",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
SHOW_DOWN <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
               "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE")
COL <- c(control = "#0072B2", caffeine = "#D55E00")                  # Okabe–Ito, CVD-safe
COL_DIR <- c(Up = "#D55E00", Down = "#0072B2", NS = "grey75")

# ---------------------------------------------------------------- inputs
s2c <- fread("metadata/samples.tsv")
s2c <- data.frame(sample = s2c$sample,
                  condition = factor(s2c$condition, levels = c("control", "caffeine")),
                  path = file.path("quant", s2c$sample), stringsAsFactors = FALSE)
s2c <- s2c[order(s2c$condition, s2c$sample), ]       # vehicle first in all figures
stopifnot(all(file.exists(file.path(s2c$path, "abundance.h5"))))
t2g <- fread("ref/t2g.tsv", data.table = FALSE)
genes <- unique(t2g[, c("ens_gene", "ext_gene", "gene_biotype")])
genes <- genes[!duplicated(genes$ens_gene), ]

# ---------------------------------------------------------------- QC (kallisto run_info.json)
qc <- rbindlist(lapply(s2c$sample, function(s) {
  j <- fromJSON(file.path("quant", s, "run_info.json"))
  data.table(sample = s, reads_processed = j$n_processed, reads_pseudoaligned = j$n_pseudoaligned,
             pct_pseudoaligned = j$p_pseudoaligned, bootstraps = j$n_bootstraps,
             kallisto_version = j$kallisto_version, call = j$call)
}))
fwrite(qc, "results/TableS1_kallisto_QC.csv")

# ---------------------------------------------------------------- sleuth (gene level)
so <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = "ens_gene", gene_mode = TRUE,
                  extra_bootstrap_summary = FALSE, read_bootstrap_tpm = FALSE,
                  num_cores = 1)
so <- sleuth_fit(so, ~condition, "full")
so <- sleuth_fit(so, ~1, "reduced")
so <- sleuth_lrt(so, "reduced", "full")
so <- sleuth_wt(so, "conditioncaffeine", "full")

lrt  <- as.data.table(sleuth_results(so, "reduced:full", "lrt", show_all = FALSE))
wald <- as.data.table(sleuth_results(so, "conditioncaffeine", "wt", show_all = FALSE))
de <- merge(wald[, .(gene_id = target_id, beta = b, se_b, mean_obs, pval_wald = pval, qval_wald = qval)],
            lrt[, .(gene_id = target_id, pval_lrt = pval, qval_lrt = qval)], by = "gene_id")
de <- merge(de, setnames(as.data.table(genes), c("gene_id", "symbol", "biotype")), by = "gene_id", all.x = TRUE)
de[, log2FC := beta / log(2)]          # sleuth beta is on the natural-log scale
de[, direction := fifelse(qval_wald <= FDR_GENE & beta > 0, "Up",
                  fifelse(qval_wald <= FDR_GENE & beta < 0, "Down", "NS"))]
setorder(de, qval_wald, pval_wald)
setcolorder(de, c("gene_id", "symbol", "biotype", "beta", "log2FC", "se_b", "mean_obs",
                  "pval_wald", "qval_wald", "pval_lrt", "qval_lrt", "direction"))
fwrite(de, "results/TableS2_differential_expression_all_genes.csv")
top_tab <- function(dir) head(de[direction == dir, .(gene_id, symbol, log2FC = round(log2FC, 2),
                                                     pval_wald = signif(pval_wald, 3), qval_wald = signif(qval_wald, 3),
                                                     pval_lrt = signif(pval_lrt, 3), qval_lrt = signif(qval_lrt, 3))], 10)
fwrite(top_tab("Up"),   "results/Table1_top_upregulated.csv")
fwrite(top_tab("Down"), "results/Table2_top_downregulated.csv")

# gene-level TPM matrix (normalised, observed)
tpm_long <- as.data.table(so$obs_norm)[, .(gene_id = target_id, sample, tpm)]
tpm <- dcast(tpm_long, gene_id ~ sample, value.var = "tpm")
tpm_m <- as.matrix(tpm[, -1]); rownames(tpm_m) <- tpm$gene_id
tpm_m <- tpm_m[, s2c$sample]
fwrite(cbind(gene_id = rownames(tpm_m), symbol = genes$ext_gene[match(rownames(tpm_m), genes$ens_gene)],
             as.data.table(round(tpm_m, 3))), "results/TableS4_gene_TPM.csv")

# ---------------------------------------------------------------- GSEA (preranked)
rank_dt <- de[!is.na(symbol) & !is.na(pval_wald)]
rank_dt[, score := -log10(pmax(pval_wald, 1e-300)) * sign(beta)]
rank_dt <- rank_dt[order(-abs(score))][!duplicated(symbol)]
ranks <- setNames(rank_dt$score, rank_dt$symbol)
ranks <- sort(ranks + rnorm(length(ranks), 0, 1e-10), decreasing = TRUE)   # break exact ties

msig <- as.data.table(msigdbr(species = "Mus musculus"))
msig <- msig[gs_cat == "H" | (gs_cat == "C2" & gs_subcat %in% c("CP:REACTOME", "CP:BIOCARTA", "CP:WIKIPATHWAYS"))]
msig[, collection := c(H = "Hallmark", "CP:REACTOME" = "Reactome", "CP:BIOCARTA" = "BioCarta",
                       "CP:WIKIPATHWAYS" = "WikiPathways")[fifelse(gs_cat == "H", "H", gs_subcat)]]
pathways <- split(msig$gene_symbol, msig$gs_name)
coll_of  <- unique(msig[, .(pathway = gs_name, collection)])

gsea <- as.data.table(fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500, eps = 0))
gsea <- merge(gsea, coll_of, by = "pathway")
setorder(gsea, padj, -NES)
gsea_out <- copy(gsea)[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
fwrite(gsea_out, "results/TableS3_GSEA_all_pathways.csv")

# ---------------------------------------------------------------- figure helpers
theme_pub <- function(base = 10) {   # minimum text size in all figures: 8 pt
  theme_classic(base_size = base, base_family = "sans") +
    theme(axis.text = element_text(colour = "black", size = base),
          axis.title = element_text(size = base + 1),
          axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
          legend.title = element_text(size = base), legend.text = element_text(size = base),
          strip.background = element_blank(), strip.text = element_text(face = "italic", size = base + 1),
          plot.title = element_text(size = base + 1, face = "bold"),
          plot.tag = element_text(size = 14, face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = "pdf")
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = w, height = h, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 200, bg = "white")
}
pretty_path <- function(x, width = 48) {
  x <- sub("^(HALLMARK|REACTOME|BIOCARTA|WP)_", "", x)
  x <- gsub("_", " ", tolower(x))
  x <- gsub("\\b(dna|rna|tnfa|nfkb|il6|jak|stat3|tca|atp|mtorc1|kras|ecm|uv|g2m|e2f|p53|tgf|tnf|nf|pi3k|akt|mtor|ppar|ros|ifng|ii|nocgmppkg)\\b",
            "\\U\\1", x, perl = TRUE)
  x <- sub("NOCGMPPKG", "NO/cGMP/PKG", x)
  x <- sub("^(.)", "\\U\\1", x, perl = TRUE)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 1), "…"), x)
}
lab_sample <- function(x) sub("caffeine_", "Caffeine ", sub("control_", "Control ", x))

# ---------------------------------------------------------------- Figure 1: PCA + volcano
logt <- log2(tpm_m + 1)
keep <- order(apply(logt, 1, var), decreasing = TRUE)[1:min(500, nrow(logt))]
pc <- prcomp(t(logt[keep, ]), center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pca_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], sample = rownames(pc$x),
                     condition = s2c$condition[match(rownames(pc$x), s2c$sample)])
fwrite(cbind(pca_df, var_PC1 = ve[1], var_PC2 = ve[2]), "results/PCA_coordinates.csv")
p1a <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = lab_sample(sample)), size = 3.4, show.legend = FALSE, seed = 1,
                  box.padding = 0.5, min.segment.length = 0) +
  scale_colour_manual(values = COL, labels = c("Vehicle", "Caffeine"), name = NULL) +
  labs(x = sprintf("PC1 (%.1f%% variance)", ve[1]), y = sprintf("PC2 (%.1f%% variance)", ve[2])) +
  theme_pub() + theme(legend.position = "top")

vol <- copy(de)[!is.na(pval_wald)]
vol[, nlp := -log10(pmax(pval_wald, 1e-300))]
lab_genes <- unique(c(head(vol[direction == "Up"], 6)$gene_id, head(vol[direction == "Down"], 6)$gene_id,
                      vol[symbol %in% KEY_GENES & direction != "NS"]$gene_id))
n_up <- sum(vol$direction == "Up"); n_dn <- sum(vol$direction == "Down")
p1b <- ggplot(vol, aes(log2FC, nlp)) +
  geom_point(data = vol[direction == "NS"], colour = COL_DIR["NS"], size = 0.6, alpha = 0.6) +
  geom_point(data = vol[direction != "NS"], aes(colour = direction), size = 1.6) +
  geom_text_repel(data = vol[gene_id %in% lab_genes], aes(label = symbol), size = 3.1, fontface = "italic",
                  max.overlaps = 30, seed = 1, min.segment.length = 0, segment.size = 0.3, box.padding = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = COL_DIR, name = sprintf("Caffeine vs vehicle (FDR < %.2f)", FDR_GENE),
                      labels = c(Down = sprintf("Down (n = %d)", n_dn), Up = sprintf("Up (n = %d)", n_up))) +
  guides(colour = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5)) +
  labs(x = expression(log[2]~"fold change (caffeine / vehicle)"), y = expression(-log[10]~italic(P)~"(Wald)")) +
  theme_pub() + theme(legend.position = "top", legend.title = element_text(size = 9, face = "bold"))
p1a <- p1a + scale_colour_manual(values = COL, labels = c("Vehicle", "Caffeine"), name = "Group") +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  theme(legend.title = element_text(size = 9, face = "bold"))
fig1 <- (p1a | p1b) + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A")
save_fig(fig1, "Figure1_PCA_volcano", 7.2, 4.2)

# ---------------------------------------------------------------- Figure 2: GSEA
sig <- gsea[padj <= FDR_PATH]
bar <- rbind(head(sig[NES > 0][order(-NES)], 8), head(sig[NES < 0][order(NES)], 8))
if (nrow(bar) > 0) {
  bar[, label := paste0(pretty_path(pathway, 40), " (", c(Hallmark = "H", Reactome = "RE", BioCarta = "BC", WikiPathways = "WP")[collection], ")")]
  bar[, label := make.unique(label, sep = " ")]
  bar[, label := factor(label, levels = label[order(NES)])]
  bar[, dir := fifelse(NES > 0, "Up", "Down")]
  p2a <- ggplot(bar, aes(NES, label, fill = dir)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = COL_DIR, name = "Enriched in",
                      labels = c(Down = "Vehicle", Up = "Caffeine")) +
    labs(x = "Normalized enrichment score (NES)", y = NULL,
         caption = sprintf("Top 8 pathways per direction, FDR < %.2f.\nH, Hallmark; RE, Reactome; BC, BioCarta; WP, WikiPathways.", FDR_PATH)) +
    theme_pub() + theme(legend.position = "top", legend.title = element_text(size = 9, face = "bold"),
                        plot.caption = element_text(size = 8, colour = "grey25", hjust = 0))
} else {
  p2a <- ggplot() + annotate("text", 0, 0, label = "No pathways at FDR < 0.05") + theme_void()
}
enr_plot <- function(pw) {
  r <- gsea[pathway == pw]
  ed <- plotEnrichmentData(pathways[[pw]], ranks)
  col <- if (r$NES > 0) COL_DIR[["Up"]] else COL_DIR[["Down"]]
  ttl <- paste(strwrap(paste0(pretty_path(pw, 80), " (", r$collection, ")"), width = 28), collapse = "\n")
  ggplot(ed$curve, aes(rank, ES)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_hline(yintercept = c(ed$posES, ed$negES)[which.max(abs(c(ed$posES, ed$negES)))],
               linewidth = 0.3, linetype = "dashed", colour = "grey40") +
    geom_segment(data = ed$ticks, aes(x = rank, xend = rank, y = -ed$spreadES / 16, yend = ed$spreadES / 16),
                 linewidth = 0.15, inherit.aes = FALSE) +
    geom_line(colour = col, linewidth = 0.7) +
    labs(title = ttl, subtitle = sprintf("NES = %.2f, FDR = %s", r$NES, formatC(r$padj, format = "e", digits = 1)),
         x = "Gene rank", y = "Enrichment score") +
    scale_x_continuous(breaks = scales::pretty_breaks(3), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_pub(8.5) + theme(plot.title = element_text(size = 8.5, face = "bold", lineheight = 0.95),
                           plot.subtitle = element_text(size = 8), axis.title = element_text(size = 9))
}
up3 <- intersect(SHOW_UP, sig[NES > 0]$pathway);   if (!length(up3)) up3 <- head(sig[NES > 0][order(padj, -NES)]$pathway, 3)
dn3 <- intersect(SHOW_DOWN, sig[NES < 0]$pathway); if (!length(dn3)) dn3 <- head(sig[NES < 0][order(padj, NES)]$pathway, 3)
row_of <- function(pws, tag) {
  if (!length(pws)) return(NULL)
  pl <- lapply(pws, enr_plot)
  while (length(pl) < 3) pl <- c(pl, list(plot_spacer()))
  wrap_elements(full = wrap_plots(pl, nrow = 1)) + labs(tag = tag)
}
rows <- Filter(Negate(is.null), list(row_of(up3, "B"), row_of(dn3, "C")))
fig2 <- wrap_elements(full = p2a) + labs(tag = "A")
if (length(rows)) fig2 <- wrap_plots(c(list(fig2), rows), ncol = 1, heights = c(2.1, rep(1.15, length(rows))))
fig2 <- fig2 & theme(plot.tag = element_text(size = 14, face = "bold"))
save_fig(fig2, "Figure2_GSEA", 7.2, 8.75)

# ---------------------------------------------------------------- Figure 3: heatmap + key genes
hm_ids <- c(head(de[direction == "Up"]$gene_id, 15), head(de[direction == "Down"]$gene_id, 15))
if (length(hm_ids) >= 2) {
  z <- t(scale(t(log2(tpm_m[hm_ids, , drop = FALSE] + 1))))
  z[is.nan(z)] <- 0
  hd <- as.data.table(z, keep.rownames = "gene_id")
  hd <- melt(hd, id.vars = "gene_id", variable.name = "sample", value.name = "z")
  hd[, symbol := make.unique(genes$ext_gene[match(hm_ids, genes$ens_gene)])[match(gene_id, hm_ids)]]
  hd[, symbol := factor(symbol, levels = rev(make.unique(genes$ext_gene[match(hm_ids, genes$ens_gene)])))]
  hd[, dir := factor(fifelse(de$direction[match(gene_id, de$gene_id)] == "Up", "Up in caffeine", "Down in caffeine"),
                     levels = c("Up in caffeine", "Down in caffeine"))]
  hd[, group := factor(fifelse(grepl("^caffeine", sample), "Caffeine", "Vehicle"), levels = c("Vehicle", "Caffeine"))]
  hd[, sample := factor(lab_sample(as.character(sample)), levels = lab_sample(s2c$sample))]
  hd[, z := pmax(pmin(z, 2), -2)]
  p3a <- ggplot(hd, aes(sample, symbol, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(dir ~ group, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                         limits = c(-2, 2), breaks = -2:2, name = "Row\nz-score") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_pub(9) +
    theme(axis.text.y = element_text(face = "italic", size = 8.5), axis.text.x = element_text(size = 9),
          axis.line = element_blank(), axis.ticks = element_blank(),
          strip.text.x = element_text(face = "bold", size = 10), strip.text.y = element_text(face = "bold", size = 9.5, angle = -90),
          legend.title = element_text(size = 9, face = "bold"), legend.key.height = unit(0.9, "cm"),
          panel.spacing = unit(3, "pt"))
} else p3a <- plot_spacer()

kg <- de[symbol %in% KEY_GENES]
kg <- kg[order(match(symbol, KEY_GENES))][!duplicated(symbol)]
if (nrow(kg)) {
  d3 <- as.data.table(tpm_m[kg$gene_id, , drop = FALSE], keep.rownames = "gene_id")
  d3 <- melt(d3, id.vars = "gene_id", variable.name = "sample", value.name = "tpm")
  d3 <- merge(d3, kg[, .(gene_id, symbol, qval_wald, log2FC)], by = "gene_id")
  d3[, condition := factor(s2c$condition[match(sample, s2c$sample)], levels = c("control", "caffeine"))]
  fmt_q <- function(q) sub("e-0?", "%*%10^-", formatC(q, format = "e", digits = 1))
  strip_lab <- setNames(sprintf("atop(italic('%s'), q == %s)", kg$symbol, fmt_q(kg$qval_wald)), kg$symbol)
  d3[, symbol := factor(strip_lab[symbol], levels = strip_lab)]
  summ <- d3[, .(mean = mean(tpm), se = sd(tpm) / sqrt(.N), ymax = max(tpm)), by = .(symbol, condition)]
  p3b <- ggplot(d3, aes(condition, tpm, colour = condition)) +
    geom_col(data = summ, aes(y = mean, fill = condition), colour = NA, alpha = 0.25, width = 0.6) +
    geom_errorbar(data = summ, aes(y = mean, ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.08, seed = 1), size = 2) +
    facet_wrap(~symbol, scales = "free_y", nrow = 1, labeller = label_parsed) +
    scale_colour_manual(values = COL, guide = "none") + scale_fill_manual(values = COL, guide = "none") +
    scale_x_discrete(labels = c("Vehicle", "Caffeine")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
    labs(x = NULL, y = "Expression (TPM)", caption = "Bars, group mean; error bars, SEM; points, individual animals; q, Wald-test FDR.") +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 10), strip.clip = "off",
                        plot.caption = element_text(size = 8, colour = "grey25", hjust = 0))
} else p3b <- plot_spacer()
fig3 <- (wrap_elements(full = p3a) + labs(tag = "A")) / (wrap_elements(full = p3b) + labs(tag = "B")) +
  plot_layout(heights = c(2.5, 1.1)) &
  theme(plot.tag = element_text(size = 14, face = "bold"))
save_fig(fig3, "Figure3_heatmap_genes", 7.2, 8.75)

# ---------------------------------------------------------------- summary for the manuscript
chk <- function(pat) {
  x <- gsea[grepl(pat, pathway)][order(padj)]
  if (!nrow(x)) return(sprintf("  %-55s not tested", pat))
  paste(sprintf("  %-55s NES = %6.2f  FDR = %.2e  (%s)", x$pathway, x$NES, x$padj,
                ifelse(x$NES > 0, "UP in caffeine", "DOWN in caffeine")), collapse = "\n")
}
sink("results/SUMMARY.txt")
cat("=== Reanalysis summary ===\n\n")
cat(readLines("results/run_parameters.txt", warn = FALSE), sep = "\n")
cat("\nR:", R.version.string, "\n")
for (p in c("sleuth", "fgsea", "msigdbr", "ggplot2", "pheatmap"))
  cat(sprintf("%s: %s\n", p, as.character(packageVersion(p))))
cat("MSigDB gene sets: msigdbr", as.character(packageVersion("msigdbr")),
    "(mouse orthologs of MSigDB human gene sets)\n")
cat("GSEA: fgseaMultilevel, minSize 15, maxSize 500, eps 0, seed 42;",
    "ranking = -log10(Wald p) * sign(beta); gene sets tested:", nrow(gsea), "\n\n")
cat("--- kallisto QC ---\n"); print(qc[, .(sample, reads_processed, pct_pseudoaligned)])
cat(sprintf("\n--- Differential expression (Wald q <= %.2f, primary) ---\nUp: %d   Down: %d   (genes tested: %d)\n",
            FDR_GENE, n_up, n_dn, nrow(de)))
cat("LRT q <= 0.05 (sensitivity): Up", sum(de$qval_lrt <= 0.05 & de$beta > 0, na.rm = TRUE),
    " Down", sum(de$qval_lrt <= 0.05 & de$beta < 0, na.rm = TRUE), "; min LRT q =", signif(min(de$qval_lrt, na.rm = TRUE), 3), "\n")
cat("Protein-coding among Wald-significant genes:", sum(de$direction != "NS" & de$biotype == "protein_coding", na.rm = TRUE), "\n")
cat("\n--- Key genes from the original manuscript ---\n")
print(de[symbol %in% KEY_GENES, .(symbol, beta = round(beta, 2), log2FC = round(log2FC, 2), qval_wald = signif(qval_wald, 3),
                                  pval_lrt = signif(pval_lrt, 3), qval_lrt = signif(qval_lrt, 3), direction)])
cat("\n--- Direction check for pathways reported in the original manuscript ---\n")
cat(chk("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|BETA_?OXIDATION|COLLAGEN|EPITHELIAL_MESENCHYMAL"), "\n")
cat(sprintf("\n--- Pathways at FDR < 0.05: %d up, %d down; FDR < 0.25: %d ---\n",
            sum(sig$NES > 0), sum(sig$NES < 0), sum(gsea$padj < 0.25)))
sink()
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
saveRDS(so, "results/sleuth_object.rds")
cat(readLines("results/SUMMARY.txt"), sep = "\n")
