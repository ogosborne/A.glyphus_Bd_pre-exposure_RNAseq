# general
library(tidyverse)
library(stringi)
# DEGs
library(DESeq2)
library(tximport)
# GO
library(GOSim)
library(topGO)
# expression-microbiome correlation
library(Hmisc) 
library(phyloseq) 
# plotting
library(ggpubr)
library(rstatix)
library(gtools)
library(gplots) 
library(GetoptLong) 
library(gridExtra)
library(grid)
library(circlize)
library(cowplot)


##### 1. PREP DATA
# load metadata
md <- read.csv("data/A.glyphus.metadata.csv", row.names = 1)
# load counts
txi <- readRDS("quant_results/txi.RDS")
# load tx2gene
tx2gene <- read.table("tr_ass_results/trinity/Aglyphus.Trinity.fasta.gene_trans_map", sep = "\t")
# switch columns
tx2gene <- tx2gene[,c(2,1)]
colnames(tx2gene) <- c("transcript", "gene")
# load protein names
prots <- read.table("tr_ass_results/transdecoder/Aglyphus.Trinity.fasta.transdecoder.gff3", sep = "\t")[,1] %>%
  stri_remove_empty_na() %>%
  unique()
# load emapper
em <- read.csv("tr_ass_results/emapper/Aglyphus.emapper.annotations", comment.char = "#", sep = "\t")
colnames(em) <- c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
em$query <- gsub("\\..*", "", em$query)
em$gene <- tx2gene[match(em$query, tx2gene[,"transcript"]),"gene"]
# filt tx2gene to keep only transcripts with proteins
tx2gene <- tx2gene[tx2gene$transcript %in% prots,]
# save gene info
dir.create("results/DESeq", recursive = T)
save(list = c("em", "tx2gene"), file = "results/DESeq/gene.info.RData")
# summarise txi to gene level (also removes transcripts filtered from tx2gene)
txi_gene <- tximport::summarizeToGene(txi, tx2gene = tx2gene)

##### 2. DIFFERENTIAL EXPRESSION
# make deseq objects, one with all samples and one for each tissue
# fun to subset txi and make dds
subset_txi2dds <- function(txi, md, col, grp, des){
  my.txi <- txi
  sel <- which(md[,col] %in% grp)
  my.txi$abundance <- txi$abundance[,sel]
  my.txi$counts <- txi$counts[,sel]
  my.txi$length <- txi$length[,sel]
  my.md <- md[sel,]
  my.dds <- DESeqDataSetFromTximport(txi = my.txi,
                                     colData = my.md,
                                     design = des)
  my.dds
}
# all data
dds <- list()
dds$all <- DESeqDataSetFromTximport(txi = txi_gene,
                                    colData = md,
                                    design = ~ Trt_TP)
# spleen data
dds$spleen <- subset_txi2dds(txi = txi_gene,
                             md = md,
                             col = "Tissue",
                             grp = "Spleen",
                             des = as.formula(~ Treatment))
# TP1 and TP2 Toe data
dds$toe <- subset_txi2dds(txi = txi_gene,
                          md = md,
                          col = "TimePoint",
                          grp = c("TP1", "TP2"),
                          des = as.formula(~ Sex + Trt_TP))

# TP1 and TP2 Toe data
#dds$toe2 <- subset_txi2dds(txi = txi_gene,
#                          md = md,
#                          col = "TimePoint",
#                          grp = c("TP1", "TP2"),
#                          des = as.formula(~ Sex + TimePoint + Treatment))
# filter out rarely expressed genes
filt_dds_lowexp <- function(dds, min.samp, min.reads = 10){
  keep <- rowSums(counts(dds) >= min.reads) >= min.samp
  my.dds <- dds[keep,]
  my.dds
}
dds$all <- filt_dds_lowexp(dds = dds$all, min.samp = 3)
dds$spleen <- filt_dds_lowexp(dds = dds$spleen, min.samp = 3)
dds$toe <- filt_dds_lowexp(dds = dds$toe, min.samp = 10)
# variance stabilising transformation
vst <- lapply(dds, vst, blind = TRUE)
save(vst, file = "results/DESeq/vst.RData")

## PCAs
dir.create("results/Ordination/")
# cols
tissue_cols <- c(Spleen = "lightcoral", Toe = "gold")
preexp_cols <- c(Naïve =  "#732f30", `Pre-exposed` = "#b38711",  Untreated = "#5b859e")
# tp shapes
tpshapes <- function(md){
  shp <- ifelse(md$Infected == "Infected", "t", "c")
  fil <- ifelse(md$TimePoint == "TP1", "o", "f")
  sf <- paste0(shp, fil)
  shapes <- c(to = 24, tf = 17, co = 21, cf = 19)
  my.shapes <- shapes[sf]
  my.shapes
}
infshapes <- function(md){
  shp <- ifelse(md$Infected == "Infected", "t", "c")
  shapes <- c(t = 24, c = 21)
  my.shapes <- shapes[shp]
  my.shapes
}
# Tissue
PCAdat <- plotPCA(vst$all,  intgroup = c("Tissue"), returnData = T) 
pdf("results/Ordination/tissue_PCA.pdf", width = 7, height = 4)
layout(matrix(c(1,2,2,2,2,3,3), ncol = 7))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     bg = tissue_cols[as.character(PCAdat$Tissue)], 
     pch = 21, 
     cex = 2,
     lwd = 0.5,
     cex.lab = 1.4,
     xlab = "PC1",
     ylab = "PC2",
     main = "",
)
plot.new()
legend("left", legend = names(tissue_cols), pt.bg = tissue_cols,  pt.cex = 3, pt.lwd = 0.5, pch = 21, bty = "n", cex = 1.6, title = expression(bold("Tissue")), title.col = "black", text.col = "gray30", title.adj = 0)
dev.off()

## Toe sample multipanel treatment + Bd
pdf("results/Ordination/treatment+Bd_PCA.pdf", width = 7, height = 8)
layout(matrix(c(1,2,2,2,2,3,3,3,
                4,5,5,5,5,6,6,6), 
              ncol = 8, byrow = T))
PCAdat <- plotPCA(vst$toe, intgroup = c("PreExposure", "Infected", "TimePoint"), returnData = T) 
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     col =  preexp_cols[as.character(PCAdat$PreExposure)], 
     pch =  tpshapes(PCAdat), 
     cex = 2,
     xlab = "PC1",
     ylab = "PC2",
     cex.lab = 1.4,
     main = "")
plot.new()
legend("topleft", 
       legend = names(preexp_cols),
       pt.cex = 3,
       col = "white",
       cex = 1.6,
       pt.bg = preexp_cols, 
       pch = 22,
       bty = "n", title = expression(bold("Treatment")),
       title.col = "black", text.col = "gray30", title.adj = 0)
legend("left", 
       legend = c("yes", "no"), 
       pt.cex = 2,
       cex = 1.6,
       xjust = 0,
       pch = c(24, 21), 
       bty = "n", title = expression(bold("Bd")),
       title.col = "black", text.col = "gray30", title.adj = 0)
legend("bottomleft", 
       legend = c("Day 14","Day 29"),
       pt.cex = 3,
       cex = 1.6,
       pch = 22,
       pt.bg = c("white", "black"),
       bty = "n", title = expression(bold("Day")),
       title.col = "black", text.col = "gray30", title.adj = 0)
# Bd load
PCAdat <- plotPCA(vst$toe, intgroup = c("logBdload","PreExposure", "Infected", "TimePoint"), returnData = T) 
Bd_cols <- colorRampPalette(c("blue4", "orangered1", "lightyellow"))
seq <- seq(from = 0,to = max(PCAdat$logBdload, na.rm = T),by = 0.1)
cols <- Bd_cols(length(seq))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     bg = cols[findInterval(PCAdat$logBdload, seq)], 
     pch =  infshapes(PCAdat),
     cex = 2,
     cex.lab = 1.4,
     xlab = "PC1",
     ylab = "PC2")
# col legend
legend_image <- as.raster(matrix(Bd_cols(20), ncol=20))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
text(x = c(0, 1.8), y = 0.75, labels = c(0,round(max(PCAdat$logBdload),2)))
rasterImage(legend_image, 0, 0.85, 1.8, 0.9)
rect(0, 0.85, 1.8, 0.9)
legend("topleft", cex = 1.6, title = expression(bold("log(Bd load)")),
       title.col = "black", legend = "", bg = rgb(0,0,0,0), bty = "n", title.adj = 0)
dev.off()

# Itraconazole
PCAdat <- plotPCA(vst$toe, intgroup = c("PreExposure", "Infected", "TimePoint", "Treatment"), returnData = T) 
PCAdat <- PCAdat[which(PCAdat$Treatment %in% c("UTA", "UVC")),]
pdf("results/Ordination/itra_PCA.pdf", width = 7, height = 4)
layout(matrix(c(1,2,2,2,2,3,3), ncol = 7))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     col =  preexp_cols[as.character(PCAdat$PreExposure)], 
     pch =  tpshapes(PCAdat), 
     cex = 2,
     xlab = "PC1",
     ylab = "PC2",
     cex.lab = 1.4,
     main = "")
plot.new()
legend("topleft", 
       legend = c("Treated", "Untreated"),
       pt.cex = 3,
       col = "white",
       cex = 1.6,
       pt.bg = preexp_cols[c(1,3)], 
       pch = 21,
       bty = "n", title = expression(bold("Itraconazole")),
       title.col = "black", text.col = "gray30", title.adj = 0)
legend("bottomleft", 
       legend = c("Day 14","Day 29"),
       pt.cex = 3,
       cex = 1.6,
       pch = 21,
       pt.bg = c("white", "black"),
       bty = "n", title = expression(bold("Day")),
       title.col = "black", text.col = "gray30", title.adj = 0)
dev.off()

# run DEseq
dds$toe <- DESeq(dds$toe)
save(dds, file = "results/DESeq/dds.RData")

# extract results
cont <- read.csv("data/contrasts.csv")
res <- list()
for(i in 1:nrow(cont)){
  gr1 <- paste(cont[i, "group1"], cont[i, "timepoint"], sep = "_")
  gr2 <- paste(cont[i, "group2"], cont[i, "timepoint"], sep = "_")
  nam <- cont[i, "test_name"]
  res[[nam]] <- results(dds$toe, contrast = c("Trt_TP", gr1, gr2))
}
# add eggnog-mapper results to results
add_em <- function(r, em){
  cols <- c("Description", "Preferred_name", "GOs")
  for(col in cols){
    r[,col] <- em[match(rownames(r), em[,"gene"]),col]
  }
  r
}
res <- lapply(res, function(x) add_em(r = x, em = em))

# combine results
res_all <- res 
# add gene column
for(i in names(res_all)){
  res_all[[i]]$gene <- rownames(res_all[[i]])
} 
# add contrast name to column names
for(i in names(res_all)){
  colnames(res_all[[i]])[c(2,5,6)] <- paste(i, colnames(res_all[[i]])[c(2,5,6)], sep = ".")
} 
# select columns and convert to data frame
for(i in names(res_all)){
  res_all[[i]] <- res_all[[i]][,c(10,7,8,9,2,5,6)] %>% as.data.frame()
} 
# merge all
res_all <- Reduce(function(x,y) merge(x, y, all = T),  res_all)

## get result intersections
mainres <- c("UV_TP1.vs.UVC_TP1","UV_TP2.vs.UVC_TP2","V_TP1.vs.VC_TP1","V_TP2.vs.VC_TP2")
int <- data.frame(matrix(NA, ncol = 4, nrow = nrow(res_all), dimnames = list(NULL, mainres)))
# over-expressed
for(i in mainres){
  int[,i] <- ifelse((res_all[,paste0(i,".padj")] >= 0.05 | is.na(res_all[,paste0(i,".padj")])) | res_all[,paste0(i,".log2FoldChange")] < 0, 
                    yes = 0, no = 1)
}
res_all$sigOvrexp_comb <- apply(int, 1, paste0, collapse="")
# under-expressed
for(i in mainres){
  int[,i] <- ifelse((res_all[,paste0(i,".padj")] >= 0.05 | is.na(res_all[,paste0(i,".padj")])) | res_all[,paste0(i,".log2FoldChange")] > 0, 
                    yes = 0, no = 1)
}
res_all$sigUndexp_comb <- apply(int, 1, paste0, collapse="")
rm(int)
# Pre-exposure specific
res_preExposure_specific <- res_all[which(res_all$sigOvrexp_comb == "0011" | res_all$sigUndexp_comb == "0011"),]
res_preExposure_specific <- res_preExposure_specific[order(apply(res_preExposure_specific[,c("V_TP1.vs.VC_TP1.padj", "V_TP2.vs.VC_TP2.padj")],1, sum), decreasing = F),]
# Naive specific
res_naive_specific <- res_all[which(res_all$sigOvrexp_comb == "1100" | res_all$sigUndexp_comb == "1100"),]
res_naive_specific <-res_naive_specific[order(apply(res_naive_specific[,c("UV_TP1.vs.UVC_TP1.padj", "UV_TP2.vs.UVC_TP2.padj")],1, sum), decreasing = F),]
# Bd specific
res_BdExp2_specific <- res_all[which(res_all$sigOvrexp_comb == "1111" | res_all$sigUndexp_comb == "1111"),]
res_BdExp2_specific <-res_BdExp2_specific[order(apply(res_BdExp2_specific[,c("UV_TP1.vs.UVC_TP1.padj","UV_TP2.vs.UVC_TP2.padj","V_TP1.vs.VC_TP1.padj","V_TP2.vs.VC_TP2.padj")],1, sum), decreasing = F),]
# write to csv
write.csv(res_preExposure_specific, file = "results/DESeq/DESeqres_preExposure_specific.csv", row.names = F)
write.csv(res_naive_specific, file = "results/DESeq/DESeqres_naive_specific.csv", row.names = F)
write.csv(res_BdExp2_specific, file = "results/DESeq/DESeqres_BdExposure2_specific.csv", row.names = F)
write.csv(res_all, file = "results/DESeq/DESeqres_all.csv", row.names = F)
# save
save(res_all, file = "results/DESeq/deseq.res_comb.RData")

##### 3. GO ENRICHMENT
ont="BP"
alg="elim"
GOanc <- GOSim::getAncestors()
# sig lists
siglists <- list()
siglists$preExp <- res_preExposure_specific$gene
siglists$noPreExp <- res_naive_specific$gene
siglists$Bd <- res_BdExp2_specific$gene
# Gene-GO term mapping
geneID2GO <- as.list(res_all$GOs)
names(geneID2GO) <- res_all$gene
geneID2GO <- lapply(geneID2GO, function(x) str_split(x, pattern = ",", simplify = F) %>% unlist())
# run for each set of DEGs
GOEnDat <- list()
GOEnTes <- list()
GOEnRes <- list()
for(i in names(siglists)){
  # list of sig genes
  my.siglist <- siglists[[i]]
  # gene list
  geneList <- ifelse(res_all$gene %in% my.siglist, yes = 1, no = 0)
  names(geneList) <- res_all$gene
  # GO datasets
  GOdata <- new("topGOdata",
                ontology = ont,
                geneSelectionFun = function(x) x > 0,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO,
                nodeSize = 10)
  GOEnDat[[i]] <- GOdata
  # run tests
  GOtests <- runTest(GOdata, algorithm = alg, statistic = "fisher")
  GOEnTes[[i]] <- GOtests
  # get results
  goRes <- GenTable(GOdata, 
                    P = GOtests,
                    topNodes = length(which(GOtests@score < 0.05)),
                    numChar = 1000)
  # add is immune
  goRes$Immune.related <- unlist(lapply(as.list(goRes$GO.ID), function(x) "GO:0002376" %in% GOanc[[x]]))
  goRes$Adaptive.Immune <- unlist(lapply(as.list(goRes$GO.ID), function(x) "GO:0002250" %in% GOanc[[x]]))
  goRes$Innate.Immune <- unlist(lapply(as.list(goRes$GO.ID), function(x) "GO:0045087" %in% GOanc[[x]]))
  # add to output
  GOEnRes[[i]] <- goRes
}
dir.create("results/GO")
save(list = c("GOEnDat","GOEnTes","GOEnRes"), file = "results/GO/GOenr.RData")
write.csv(GOEnRes$preExp[,1:6], file = "results/GO/GOEnrich_preExposure_specific.csv", row.names = F)
write.csv(GOEnRes$noPreExp[,1:6], file = "results/GO/GOEnrich_naive_specific.csv", row.names = F)
write.csv(GOEnRes$Bd[,1:6], file = "results/GO/GOEnrich_Bd_specific.csv", row.names = F)
# Remove unnecessary cols and combine all GO enrichment results
all_GO <- GO_en[["results"]]
for(i in names(all_GO)) colnames(all_GO[[i]])[6] <- paste(i, colnames(all_GO[[i]])[6], sep = ".")
for(i in names(all_GO)) all_GO[[i]] <- all_GO[[i]][,c(1,2,6,7)]
all_GO <- lapply(all_GO, as.data.frame)
all_GO <- Reduce(function(x,y) merge(x, y, all = T),  all_GO)
all_GO <- all_GO %>% replace(is.na(.), "n.s")
# put immune related terms first
all_GO <- rbind(all_GO[all_GO$Immune.related,],all_GO[!(all_GO$Immune.related),])
# save
write.csv(all_GO, "results/GO/GOenrichment.csv", row.names = F)

##### 6. CURATED GENES - MICROBIOME CORRELATION 
# load MB data
# code to create microbiome input files is in a separate repo, see publication for details 
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
# get curated genes
curated_immune <- read.csv("data/manually_curated_immune_genes.csv")
samp_select <- gsub("_","-",gsub("s","",sample_names(ps_16S)))
dat <- assay(vst$toe)[curated_immune$gene, samp_select] %>% as.data.frame()
rownames(dat) <- make.names(curated_immune$Preferred_name_short, unique = TRUE)
colnames(dat) <- sample_names(ps_16S)
dat <- t(dat)
# merge to genus
ps_16S_glom <- tax_glom(ps_16S, taxrank = "Genus")
ps_ITS_glom <- tax_glom(ps_ITS, taxrank = "Genus")
# filter taxa in less than 50% samples
ps_16S_glom_f <- filter_taxa(ps_16S_glom, function (x) {sum(x > 0) >= nsamples(ps_16S_glom)/2}, prune=TRUE)
ps_ITS_glom_f <- filter_taxa(ps_ITS_glom, function (x) {sum(x > 0) >= nsamples(ps_ITS_glom)/6}, prune=TRUE)
ntaxa(ps_16S_glom_f) ; ntaxa(ps_ITS_glom_f)
# get matrices of proportional abundance
gtab_16S <- t(as(otu_table(ps_16S_glom_f), "matrix"))
colnames(gtab_16S) <- tax_table(ps_16S_glom_f)[,"Genus"]
gtab_16S <- proportions(gtab_16S, margin = 2)
gtab_16S <- gtab_16S[row.names(dat),]
gtab_ITS <- t(as(otu_table(ps_ITS_glom_f), "matrix"))
colnames(gtab_ITS) <- tax_table(ps_ITS_glom_f)[,"Genus"]
gtab_ITS <- proportions(gtab_ITS, margin = 2)
gtab_ITS <- gtab_ITS[row.names(dat),]
# combine gene expression and bacterial genera
all_dat <- cbind(dat, gtab_16S)
all_dat <- cbind(all_dat, gtab_ITS) 
all_dat <- cbind(all_dat, metadata[rownames(gtab_16S),"logBdload"]) %>% as.matrix()
get_type_cur <- curated_immune$Category
names(get_type_cur) <- make.names(curated_immune$Preferred_name_short, unique = TRUE)
nodetype <- c(unname(get_type_cur[colnames(dat)]), rep("Bacteria", ncol(gtab_16S)), rep("Fungi", ncol(gtab_ITS)), "Bd")
names(nodetype) <- colnames(all_dat)
# get correlations
all_cors <- rcorr(all_dat, type = "spearman")
# adjust P vals
all_cors$Q <- p.adjust(all_cors$P[,], method = "fdr") %>% 
  matrix(ncol = ncol(all_cors$P), dimnames = list(rownames(all_cors$P), colnames(all_cors$P)))
