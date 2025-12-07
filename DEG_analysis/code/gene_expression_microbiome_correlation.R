# Microbiome
library(phyloseq) # phyloseq
# Expression
library(DESeq2)
# general
library(tidyverse)
# sparse CCA
library(compositions) # clr correction
library(PMA) # cca
# gene-microbe correlations
library(Hmisc) # rcorr
library(cocor) # z test
library(PermCor)
# plotting
library(ggplot2)
library(patchwork)
# load functions
source("DEG_analysis/code/gene_expression_microbiome_correlation_funcs.R")
# outdir
outdir="results/CCA_Results/"
dir.create(outdir)

##### 1. PREP DATA
# load MB data
# ASV table
featfile_16S <- "data/mb/featureTab_16S.csv"
featureTab_16S <- as.matrix(read.csv(featfile_16S, header = T, row.names = 1))
featureTab_16S <- otu_table(featureTab_16S, taxa_are_rows = TRUE)
featfile_ITS <- "data/mb/featureTab_ITS.csv"
featureTab_ITS <- as.matrix(read.csv(featfile_ITS, header = T, row.names = 1))
featureTab_ITS <- otu_table(featureTab_ITS, taxa_are_rows = TRUE)
# taxonomy
taxofile_16S <- "data/mb/taxonomy_16S.csv"
taxonomy_16S <- as.matrix(read.csv(taxofile_16S, row.names = 1))
taxonomy_16S <- tax_table(taxonomy_16S)
taxofile_ITS <- "data/mb/taxonomy_ITS.csv"
taxonomy_ITS <- as.matrix(read.csv(taxofile_ITS, row.names = 1))
taxonomy_ITS <- tax_table(taxonomy_ITS)
# metadata
metafile <- "data/mb/microbiome.MD.csv"
metadata <- read.csv(metafile, header = T, row.names = 1)
metadata <- sample_data(metadata)
# phyloseq
ps_16S <- merge_phyloseq(featureTab_16S, taxonomy_16S, metadata) 
ps_ITS <- merge_phyloseq(featureTab_ITS, taxonomy_ITS, metadata) 
# split timepoints
ps_16S_TP1 <- subset_samples(ps_16S, TimePoint == "TP1")
ps_16S_TP2 <- subset_samples(ps_16S, TimePoint == "TP2")
ps_ITS_TP1 <- subset_samples(ps_ITS, TimePoint == "TP1")
ps_ITS_TP2 <- subset_samples(ps_ITS, TimePoint == "TP2")
# filter out taxa in less than 20% samples in either TP
# 16S
ps_16S_TP1_f <- filter_taxa(ps_16S_TP1, function (x) {sum(x > 0) >= nsamples(ps_16S_TP1)/5}, prune=TRUE)
ps_16S_TP2_f <- filter_taxa(ps_16S_TP2, function (x) {sum(x > 0) >= nsamples(ps_16S_TP2)/5}, prune=TRUE)
tax2keep_16S <- intersect(taxa_names(ps_16S_TP1_f), taxa_names(ps_16S_TP2_f))
ps_16S_TP1_f <- prune_taxa(tax2keep_16S, ps_16S_TP1_f) 
ps_16S_TP2_f <- prune_taxa(tax2keep_16S, ps_16S_TP2_f) 
# ITS
ps_ITS_TP1_f <- filter_taxa(ps_ITS_TP1, function (x) {sum(x > 0) >= nsamples(ps_ITS_TP1)/10}, prune=TRUE)
ps_ITS_TP2_f <- filter_taxa(ps_ITS_TP2, function (x) {sum(x > 0) >= nsamples(ps_ITS_TP2)/10}, prune=TRUE)
tax2keep_ITS <- intersect(taxa_names(ps_ITS_TP1_f), taxa_names(ps_ITS_TP2_f))
ps_ITS_TP1_f <- prune_taxa(tax2keep_ITS, ps_ITS_TP1_f) 
ps_ITS_TP2_f <- prune_taxa(tax2keep_ITS, ps_ITS_TP2_f) 
# rename ASVs to avoid merging bacterial and fungal ASVs
taxa_names(ps_16S_TP1_f) <- paste0("B_",taxa_names(ps_16S_TP1_f))
taxa_names(ps_16S_TP2_f) <- paste0("B_",taxa_names(ps_16S_TP2_f))
taxa_names(ps_ITS_TP1_f) <- paste0("F_",taxa_names(ps_ITS_TP1_f))
taxa_names(ps_ITS_TP2_f) <- paste0("F_",taxa_names(ps_ITS_TP2_f))
# merge 16S and ITS
ps_TP1 <-  merge_phyloseq(ps_16S_TP1_f, ps_ITS_TP1_f)
ps_TP2 <-  merge_phyloseq(ps_16S_TP2_f, ps_ITS_TP2_f)
# Extract ASV table
mb_TP1 <- as(otu_table(ps_TP1), "matrix") %>% t()
mb_TP2 <- as(otu_table(ps_TP2), "matrix") %>% t()
# rename samples to match expression data
rownames(mb_TP1) <- sub('.', '', rownames(mb_TP1))
rownames(mb_TP1) <- gsub("_", "-", rownames(mb_TP1))
rownames(mb_TP2) <- sub('.', '', rownames(mb_TP2))
rownames(mb_TP2) <- gsub("_", "-", rownames(mb_TP2))
# CLR transform (log-transformed relative abundance, centered/scaled)
mb_TP1_clr <- clr(mb_TP1 + 1e-6) %>% scale(., center = T, scale = T)
mb_TP2_clr <- clr(mb_TP2 + 1e-6) %>% scale(., center = T, scale = T)
# load expression data (previously VST transformed)
load("results/DESeq/vst.RData")
vst <- vst$toe
# subset samples by timepoint
vst_TP1 <- subset(vst, select=colData(vst)$TimePoint == "TP1")
vst_TP2 <- subset(vst, select=colData(vst)$TimePoint == "TP2")
# filter genes to keep those in pre-exposure specific and infection specific sets
load("results/DESeq/DESeq_all_contrasts.RData")
DEGs <- res_all[which(res_all$preExposure_specific | res_all$BdExp2_specific), "gene"]
vst_TP1 <- vst_TP1[DEGs,]
vst_TP2 <- vst_TP2[DEGs,]
# load curated immune genes for reference
immune <- read.csv("data/Leon_categorised_immune_genes_update.csv")
# subset samples which are in both microbiome + expression datasets
samps_TP1 <- intersect(rownames(mb_TP1_clr), colnames(vst_TP1))
samps_TP2 <- intersect(rownames(mb_TP2_clr), colnames(vst_TP2))
mb_TP1_clr <- mb_TP1_clr[samps_TP1,]
mb_TP2_clr <- mb_TP2_clr[samps_TP2,]
vst_TP1 <- vst_TP1[,samps_TP1]
vst_TP2 <- vst_TP2[,samps_TP2]
# get scaled expression matrix
expr_TP1 <- assay(vst_TP1) %>% as.matrix() %>% t() %>% scale(., center = T, scale = T)
expr_TP2 <- assay(vst_TP2) %>% as.matrix() %>% t() %>% scale(., center = T, scale = T)
# sanity check
all(rownames(mb_TP1_clr) == rownames(expr_TP1))
all(rownames(mb_TP2_clr) == rownames(expr_TP2))

##### 2. SPARSE CCA
# build data object for CCA
CCAdat <- list()
TPs <- c("TP1","TP2")
group_info <- colData(vst)
for(i in TPs){
  CCAdat[[i]][["mb"]]  <- get(paste0("mb_", i, "_clr"))
  CCAdat[[i]][["ex"]]  <- get(paste0("expr_", i))
  CCAdat[[i]][["trt"]] <- group_info[rownames(CCAdat[[i]][["ex"]]), "Treatment"]
  CCAdat[[i]][["logBdload"]] <- group_info[rownames(CCAdat[[i]][["ex"]]), "logBdload"]
}
# Residualise every column of CCAdat[[TP]]$ex and CCAdat[[TP]]$mb by logBdload. Use the residuals for CCA & correlations.
for(i in TPs){
  CCAdat[[i]][["mb_res"]]  <- apply(CCAdat[[i]][["mb"]], 2, function(y) resid(lm(y ~ CCAdat[[i]][["logBdload"]])))
  CCAdat[[i]][["ex_res"]]  <- apply(CCAdat[[i]][["ex"]], 2, function(y) resid(lm(y ~ CCAdat[[i]][["logBdload"]])))
}
# get tuning parameters
perm_res <- list()
for(i in TPs){
  perm_res[[i]] <- CCA.permute(CCAdat[[i]]$ex_res, CCAdat[[i]]$mb_res, typex="standard", typez="standard", nperms=1000)
}
# Run sparse CCA with bootstrap. Record selection frequency for each feature (i.e. the proportion of bootstraps in which the feature has a non-zero weight)
bs_cca_res_results <- list()
seeds <- c(TP1 = 769769, TP2 = 12189)
for(tp in TPs){
  bs_cca_res_results[[tp]] <- bootstrap_cca(ex = CCAdat[[tp]]$ex_res,
                                            mb = CCAdat[[tp]]$mb_res,
                                            perm_obj = perm_res[[tp]],
                                            K = 1, nboot = 1000,
                                            seed = seeds[[tp]])
}
# output table 
combined_bscca_table <- combine_bscca(bs_cca_res_results)
write.csv(combined_bscca_table, paste0(outdir,"CCA_results.csv"), row.names = FALSE)
# Plot canonical variates (with mean canonical weights across bootstrap reps)
aes <- read.csv("data/aesthetics.csv", row.names = 1, header = T)
score_plots <- list()
for(tp in TPs){
  # Mean canonical weights
  u_mean <- bs_cca_res_results[[tp]]$ex$mean_signed_weight
  names(u_mean) <- bs_cca_res_results[[tp]]$ex$feature
  v_mean <- bs_cca_res_results[[tp]]$mb$mean_signed_weight
  names(v_mean) <- bs_cca_res_results[[tp]]$mb$feature
  # Canonical variates
  u_scores <- as.matrix(CCAdat[[tp]]$ex_res) %*% u_mean
  v_scores <- as.matrix(CCAdat[[tp]]$mb_res) %*% v_mean
  # Prepare plotting data
  df <- data.frame(
    Host = u_scores[,1],
    Microbe = v_scores[,1],
    Treatment = CCAdat[[tp]]$trt
  )
  # cor test
  ct <- cor.test(u_scores, v_scores, method = "pearson")
  print(ct)
  # Build label for annotation
  r_txt <- sprintf("r = %.2f", ct$estimate)
  p_txt   <- ifelse(ct$p.value < 0.001,
                    "P < 0.001",
                    sprintf("P = %.3f", ct$p.value))
  label_txt <- paste(r_txt, p_txt, sep = "\n")
  # get aesthetics
  aes_tp <- aes[grep(tp, rownames(aes)), ]
  rownames(aes_tp) <- gsub(paste0("_", tp), "", rownames(aes_tp))
  # plot
  score_plots[[tp]] <- ggplot(df, aes(x = Host, y = Microbe)) +
    geom_point(size = 3, aes(color = Treatment, shape = Treatment)) +
    labs(title = tp, x = "Mean Canonical Variate (Host Expression)", y = "Mean Canonical Variate (Microbiome)") +
    scale_color_manual(values = setNames(aes_tp$col, rownames(aes_tp))) +
    scale_shape_manual(values = setNames(aes_tp$shp, rownames(aes_tp))) +
    geom_smooth(method='lm', formula = y ~ x, color = "black") +
    annotate("text",
             x = Inf, y = -Inf,
             hjust = 1.1, vjust = -0.5,
             label = label_txt, size = 5) +
    theme_minimal()
}
combined_plot <- score_plots$TP1 / score_plots$TP2
pdf(paste0(outdir,"canonical_variate_plots.pdf"), width = 6.75, height = 8.2)
combined_plot
dev.off()

##### 3. ALL-BY-ALL CORRELATIONS OF CCA-SELECTED FEATURES
# select features for further analysis as the features that have a bootstrap selection frequency > 0.4 in either timepoint 
thresh = 0.4
uni_ex <- union(bs_cca_res_results$TP1$ex[which(bs_cca_res_results$TP1$ex$selection_freq > thresh),"feature"],
                    bs_cca_res_results$TP2$ex[which(bs_cca_res_results$TP2$ex$selection_freq > thresh),"feature"])
uni_mb <- union(bs_cca_res_results$TP1$mb[which(bs_cca_res_results$TP1$mb$selection_freq > thresh),"feature"],
                    bs_cca_res_results$TP2$mb[which(bs_cca_res_results$TP2$mb$selection_freq > thresh),"feature"])
uni.dat <- list(TP1 = list(), TP2 = list())
for(tp in TPs){
  uni.dat[[tp]]$ex <- CCAdat[[tp]]$ex_res[,uni_ex]
  uni.dat[[tp]]$mb <- CCAdat[[tp]]$mb_res[,uni_mb]
}
# calculate all-by-all correlations
cors_uni_all <- allcors_ex_mb(dat = uni.dat, only_mb_ex = F)
# convert to data frame
corsdf_uni_all <- cors2df(cors = cors_uni_all, gene_info = res_all, taxonomy = tax_table(ps_TP1), only_mb_ex = F, symmetric = T, inc.asv = T, sig_filter = "either")
# keep only cors which are significant at both timepoints for network plot and group-wise correlation analysis
corsdf_uni_all_bothTP_allCor <- filter_by_significance(corsdf_uni_all, sig_filter = "both", alpha = 0.05, min_abs_cor = 0)
write.csv(corsdf_uni_all_bothTP_allCor, paste0(outdir,"cors_bothTPs.csv"), row.names = FALSE)
# Make input for gephi for network plot
# Edge table (gene–microbe only, with mean and sign) ---
gephi_edges <- corsdf_uni_all_bothTP_allCor %>%
  filter(cor_type == "gene-microbe") %>%
  rowwise() %>%
  mutate(
    mean_r = mean(c(TP1_r, TP2_r), na.rm = TRUE),
    sign = case_when(
      mean_r > 0 ~ 1,
      mean_r < 0 ~ -1,
      TRUE ~ 0
    )
  ) %>%
  ungroup() %>%
  select(source = feature1, target = feature2, mean_r, sign)

# Count the number of GM edges per feature 
gm_count1 <- gephi_edges  %>%
  group_by(source) %>%
  summarise(n_gm_edges = n(), .groups = "drop") %>%
  rename(id = source)

gm_count2 <- gephi_edges %>%
  group_by(target) %>%
  summarise(n_gm_edges = n(), .groups = "drop") %>%
  rename(id = target)

gm_counts <- bind_rows(gm_count1, gm_count2) %>%
  group_by(id) %>%
  summarise(n_gm_edges = sum(n_gm_edges), .groups = "drop")

# Build node table
nodes <- corsdf_uni_all_bothTP_allCor %>%
  select(id = feature1, type = feature1_type, label = feature1_name) %>%
  bind_rows(
    corsdf_uni_all_bothTP_allCor %>%
      select(id = feature2, type = feature2_type, label = feature2_name)
  ) %>%
  distinct(id, .keep_all = TRUE) %>%
  left_join(gm_counts, by = "id") %>%
  mutate(n_gm_edges = ifelse(is.na(n_gm_edges), 0, n_gm_edges))
# Keep only nodes with >=1 GM edge
nodes_filt <- nodes %>%
  filter(n_gm_edges >= 1)
# keep only edges between retained nodes
edges_filt <- gephi_edges %>%
  filter(source %in% nodes_filt$id, target %in% nodes_filt$id, abs(mean_r) > 0.1)
# remove orphan nodes
connected_ids <- union(edges_filt$source, edges_filt$target)
nodes_filt <- nodes_filt %>%
  filter(id %in% connected_ids)
# add weight
edges_filt$Weight <- abs(edges_filt$mean_r)
# add gene type
nodes_filt$type2 <- ifelse(nodes_filt$type == "gene",
                           ifelse(nodes_filt$id %in% immune$gene,
                                  immune[match(nodes_filt$id, immune$gene), "Category2"],
                                  "non-immune"),
                           nodes_filt$type)

# write to CSV for Gephi
write.csv(nodes_filt, paste0(outdir,"gephi_nodes.csv"), row.names = FALSE)
write.csv(edges_filt, paste0(outdir,"gephi_edges.csv"), row.names = FALSE)
# plot barplot
bact_nodes <- nodes_filt[nodes_filt$type == "bacteria",]
bact_nodes <- bact_nodes[order(bact_nodes$n_gm_edges, decreasing = T),]
bact_nodes$Phylum <- c(tax_table(ps_TP1)[bact_nodes$id, "Phylum"])
bac_cols <- c("#CC99BA", "#FFBFE9", "#FFDFF4")
names(bac_cols) <- sort(unique(bact_nodes$Phylum))
par(mar = c(4,15,1,4))
barplot(rev(bact_nodes$n_gm_edges[1:20]), 
        names.arg = rev(bact_nodes$label[1:20]), 
        horiz = T, 
        col = bac_cols[rev(bact_nodes$Phylum[1:20])], 
        las = 2, 
        xlim = c(0,200),
        xlab = "N gene-microbe edges")

##### 4. COMPARE GENE-MICROBE CORRELATIONS BETWEEN TREATMENT GROUPS
# keep only gene-microbe cors which are significant at both timepoints to look at differences between groups 
corsdf_uni_all_bothTP_gmOnly <- corsdf_uni_all_bothTP_allCor %>%
  filter(cor_type == "gene-microbe") 
# calculate gene-microbe correlations among samples from each treatment group
groupcors <- compute_group_cors(CCAdat = CCAdat, feature_pairs = corsdf_uni_all_bothTP_gmOnly)
# run fisher's exact tests to test whether the signs of gene-microbe correlations are significantly associated between treatment groups
sign_test <- sign_concordance_pairwise(group_corrs = groupcors)
# run permutation-based correlation tests to test whether the values of gene-microbe correlations are significantly associated between treatment groups
set.seed(981973)
cor_of_cor_test <- permtest_group_similarity(group_corrs = groupcors, B = 1000, method = "Pearson")
# add sign test results
cor_of_cor_test <- cbind(cor_of_cor_test, sign_test[,3:ncol(sign_test)])
# run fisher's Z tests to detect individual gene-microbe correlations that are significantly different between treatment groups 
ztests <- run_pairwise_cocor_tests(group_corrs = groupcors)
# plot correlations between the values of gene-microbe correlations between each pair of treatment groups
pdf(paste0(outdir,"group_correlation_matrix.pdf"), width = 10, height = 10)
plot_group_cors_matrix(group_corrs = groupcors, perm_results = cor_of_cor_test, aes = aes)
dev.off()
# write table
write.csv(cor_of_cor_test, file = paste0(outdir,"exp-mb_cors_group_comparison.csv"), row.names = F)
