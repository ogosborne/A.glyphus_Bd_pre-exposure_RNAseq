# general
library(tidyverse) # %>%
library(stringi) # stri_remove_empty_na
# DEGs
library(DESeq2)
library(tximport)
# upset plot
library(UpSetR)
# GO
library(topGO) # topGO
# expression pattern plots
library(ggpubr)
library(ggplot2)
library(rstatix)

##### 1. PREP DATA
# load metadata
md <- read.csv("data/A.glyphus.metadata.csv", row.names = 1)
# load expression counts
txi <- readRDS("quant_results/txi.RDS")
# load gene_trans_map from trinity - used to merge transcripts to gene level
tx2gene <- read.table("tr_ass_results/trinity/Aglyphus.Trinity.fasta.gene_trans_map", sep = "\t")
# switch gene_trans_map columns
tx2gene <- tx2gene[,c(2,1)]
colnames(tx2gene) <- c("transcript", "gene")
# load protein names from transdecoder - used to filter out genes with no detected ORF
prots <- read.table("tr_ass_results/transdecoder/Aglyphus.Trinity.fasta.transdecoder.gff3", sep = "\t")[,1] %>%
  stri_remove_empty_na() %>%
  unique()
# load emapper functional annotation results
em <- read.csv("tr_ass_results/emapper/Aglyphus.emapper.annotations", comment.char = "#", sep = "\t")
colnames(em) <- c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
em$query <- gsub("\\..*", "", em$query)
em$gene <- tx2gene[match(em$query, tx2gene[,"transcript"]),"gene"]
# filter gene_trans_map to keep only transcripts with with detected ORF - this will automatically filter genes when we summarise to gene level
tx2gene <- tx2gene[tx2gene$transcript %in% prots,]
# save gene info
dir.create("results/DESeq", recursive = T)
save(list = c("em", "tx2gene"), file = "results/DESeq/gene.info.RData")
# summarise counts to gene level (also removes transcripts filtered from gene_trans_map)
txi_gene <- tximport::summarizeToGene(txi, tx2gene = tx2gene)
# make deseq objects, one with all samples and one for each tissue
# function to subset txi and make dds
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

##### 2. EXPRESSION ORDINATIONS
dir.create("results/Ordination/")
# cols
aes <- read.csv("data/aesthetics.csv", row.names = 1, header = T)
# Tissue
PCAdat <- plotPCA(vst$all,  intgroup = c("Tissue"), returnData = T) 
pdf("results/Ordination/tissue_PCA.pdf", width = 7, height = 4)
layout(matrix(c(1,2,2,2,2,3,3), ncol = 7))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     bg = aes[PCAdat$group,"col"], 
     pch = 21, 
     cex = 2,
     lwd = 0.5,
     cex.lab = 1.4,
     xlab = "PC1",
     ylab = "PC2",
     main = "",
)
plot.new()
legend("left", legend = c("Spleen","Toe"), pt.bg = aes[c("Spleen","Toe"),"col"],  pt.cex = 3, pt.lwd = 0.5, pch = 21, bty = "n", cex = 1.6, title = expression(bold("Tissue")), title.col = "black", text.col = "gray30", title.adj = 0)
dev.off()
# Toe samples: multipanel treatment + Bd
pdf("results/Ordination/treatment+Bd_PCA.pdf", width = 7, height = 8)
layout(matrix(c(1,2,2,2,2,3,3,3,
                4,5,5,5,5,6,6,6), 
              ncol = 8, byrow = T))
PCAdat <- plotPCA(vst$toe, intgroup = c("Trt_TP"), returnData = T) 
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     col =  aes[PCAdat$group,"col"], 
     pch =  aes[PCAdat$group,"shp"], 
     cex = 2,
     xlab = "PC1",
     ylab = "PC2",
     cex.lab = 1.4,
     main = "")
plot.new()
grp <- PCAdat$group %>% unique() %>% as.character() %>% sort()
# legend
legend("topleft", 
       legend = gsub("_", ": ", grp),
       pt.cex = 3,
       #col = "white",
       cex = 1.6,
       col = aes[grp,"col"], 
       pch = aes[grp,"shp"],
       bty = "n", title = expression(bold("Treatment")),
       title.col = "black", text.col = "gray30", title.adj = 0)
# Bd load
PCAdat <- plotPCA(vst$toe, intgroup = c("logBdload","Trt_TP"), returnData = T) 
Bd_cols <- colorRampPalette(c("blue4", "orangered1", "lightyellow"))
seq <- seq(from = 0,to = max(PCAdat$logBdload, na.rm = T),by = 0.1)
cols <- Bd_cols(length(seq))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     bg = cols[findInterval(PCAdat$logBdload, seq)],
     col = "black",
     pch = aes[PCAdat$Bd2,"shp"],
     cex = 2,
     cex.lab = 1.4,
     xlab = "PC1",
     ylab = "PC2")
# legend
legend_image <- as.raster(matrix(Bd_cols(20), ncol=20))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
text(x = c(0, 1.8), y = 0.75, labels = c(0,round(max(PCAdat$logBdload),2)))
rasterImage(legend_image, 0, 0.85, 1.8, 0.9)
rect(0, 0.85, 1.8, 0.9)
legend("topleft", cex = 1.6, title = expression(bold("log(Bd load)")),
       title.col = "black", legend = "", bg = rgb(0,0,0,0), bty = "n", title.adj = 0)
dev.off()
# Itraconazole
PCAdat <- plotPCA(vst$toe, intgroup = c("Trt_TP"), returnData = T) 
PCAdat <- PCAdat[which(PCAdat$Treatment %in% c("untreated", "mock-mock")),]
pdf("results/Ordination/itra_PCA.pdf", width = 7, height = 4)
layout(matrix(c(1,2,2,2,2,3,3,3), ncol = 8))
plot.new()
plot(PCAdat$PC1, PCAdat$PC2, 
     col = aes[PCAdat$group,"col"],
     pch =  aes[PCAdat$group,"shp"], 
     cex = 2,
     xlab = "PC1",
     ylab = "PC2",
     cex.lab = 1.4,
     main = "")
plot.new()
grp <- PCAdat$group %>% unique() %>% as.character() %>% sort()
# legend
legend("topleft", 
       legend = gsub("_", ": ", grp),
       pt.cex = 3,
       #col = "white",
       cex = 1.6,
       col = aes[grp,"col"], 
       pch = aes[grp,"shp"],
       bty = "n", title = expression(bold("Treatment")),
       title.col = "black", text.col = "gray30", title.adj = 0)
dev.off()
##### 3. DIFFERENTIAL EXPRESSION
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
# add eggnog-mapper annotations to results
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
# upset plot to show intersection of results between tests
tests <- grep("padj", colnames(res_all), value = T) 
uslist <- list()
for(i in tests){
  nicename <- gsub(".padj", "", i) %>% 
    gsub(".vs.", " Vs ", ., fixed = T) %>% 
    gsub(".", "-", ., fixed = T) %>% 
    gsub("_", ": ", ., fixed = T) 
  uslist[[nicename]] <- res_all[which(res_all[,i] < 0.05), "gene"]
}
# upset for all contrasts
upset(fromList(uslist), 
      nsets = 10, 
      nintersects = 1000, 
      keep.order = T, 
      sets = rev(names(uslist)),
      order.by = "freq", 
      text.scale = 1.2,
      set_size.scale_max = max(lengths(uslist))*1.2)
# upset for largest four contrasta
uslist <- uslist[1:4]
upset(fromList(uslist), 
      nsets = 10, 
      nintersects = 1000, 
      keep.order = T, 
      sets = rev(names(uslist)),
      order.by = "freq", 
      text.scale = 2,
      set_size.scale_max = max(lengths(uslist))*1.2,
      query.legend = "top",
      queries = list(list(query = intersects, params = list("Bd-Bd Vs Bd-mock: TP1", "Bd-Bd Vs Bd-mock: TP2"), color = "red", active = TRUE, query.name = "Preexposure-specific"),
                     list(query = intersects, params = list("mock-Bd Vs mock-mock: TP1", "mock-Bd Vs mock-mock: TP2", "Bd-Bd Vs Bd-mock: TP1", "Bd-Bd Vs Bd-mock: TP2"), color= "blue", active = TRUE, query.name = "Bd-specific")))
## get result intersections
mainres <- c("mock.Bd.vs.mock.mock_TP1",
             "mock.Bd.vs.mock.mock_TP2",
             "Bd.Bd.vs.Bd.mock_TP1",
             "Bd.Bd.vs.Bd.mock_TP2")
int <- data.frame(matrix(NA, ncol = 4, nrow = nrow(res_all), dimnames = list(NULL, mainres)))
# separate over- and under-expression to ensure significant DE for each timepoint is in the same direction
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
# Naive specific
res_naive_specific <- res_all[which(res_all$sigOvrexp_comb == "1100" | res_all$sigUndexp_comb == "1100"),]
res_naive_specific <-res_naive_specific[order(apply(res_naive_specific[,c("mock.Bd.vs.mock.mock_TP1.padj", "mock.Bd.vs.mock.mock_TP2.padj")],1, sum), decreasing = F),]
# Pre-exposure specific
res_preExposure_specific <- res_all[which(res_all$sigOvrexp_comb == "0011" | res_all$sigUndexp_comb == "0011"),]
res_preExposure_specific <- res_preExposure_specific[order(apply(res_preExposure_specific[,c("Bd.Bd.vs.Bd.mock_TP1.padj", "Bd.Bd.vs.Bd.mock_TP2.padj")],1, sum), decreasing = F),]
# Bd specific
res_BdExp2_specific <- res_all[which(res_all$sigOvrexp_comb == "1111" | res_all$sigUndexp_comb == "1111"),]
res_BdExp2_specific <-res_BdExp2_specific[order(apply(res_BdExp2_specific[,c("mock.Bd.vs.mock.mock_TP1.padj", "mock.Bd.vs.mock.mock_TP2.padj","Bd.Bd.vs.Bd.mock_TP1.padj", "Bd.Bd.vs.Bd.mock_TP2.padj")],1, sum), decreasing = F),]
# write to csv
write.csv(res_preExposure_specific, file = "results/DESeq/DESeqres_preExposure_specific.csv", row.names = F)
write.csv(res_naive_specific, file = "results/DESeq/DESeqres_naive_specific.csv", row.names = F)
write.csv(res_BdExp2_specific, file = "results/DESeq/DESeqres_BdExposure2_specific.csv", row.names = F)
# add to main results
res_all$preExposure_specific <- res_all$gene %in% res_preExposure_specific$gene
res_all$naive_specific <- res_all$gene %in% res_naive_specific$gene
res_all$BdExp2_specific <- res_all$gene %in% res_BdExp2_specific$gene
res_all$sigOvrexp_comb <- NULL
res_all$sigUndexp_comb <- NULL
# save combined results
write_csv(res_all, file = "results/DESeq/DESeq_all_contrasts.csv", col_names = T)
save(res_all, file = "results/DESeq/DESeq_all_contrasts.RData")

##### 3. GO ENRICHMENT
ont="BP"
alg="elim"
# intersections of interest
siglists <- list()
siglists$preExp <- res_preExposure_specific$gene
siglists$naive <- res_naive_specific$gene
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
                    topNodes = length(GOtests@score),
                    numChar = 1000)
  # adjust P values
  goRes$Q <- p.adjust(goRes$P, method = "fdr")
  goRes <- goRes[which(goRes$Q < 0.05), ]
  # add to output
  GOEnRes[[i]] <- goRes
}
dir.create("results/GO")
save(list = c("GOEnDat","GOEnTes","GOEnRes"), file = "results/GO/GOenr.RData")
write.csv(GOEnRes$preExp[,1:6], file = "results/GO/GOEnrich_preExposure_specific.csv", row.names = F)
write.csv(GOEnRes$naive[,1:6], file = "results/GO/GOEnrich_naive_specific.csv", row.names = F)
write.csv(GOEnRes$Bd[,1:6], file = "results/GO/GOEnrich_Bd_specific.csv", row.names = F)
# Remove unnecessary cols and combine all GO enrichment results
all_GO <- GOEnRes
for(i in names(all_GO)) colnames(all_GO[[i]])[6] <- paste(i, colnames(all_GO[[i]])[6], sep = ".")
for(i in names(all_GO)) colnames(all_GO[[i]])[7] <- paste(i, colnames(all_GO[[i]])[7], sep = ".")
for(i in names(all_GO)) all_GO[[i]] <- all_GO[[i]][,c(1,2,6,7)]
all_GO <- lapply(all_GO, as.data.frame)
all_GO <- Reduce(function(x,y) merge(x, y, all = T),  all_GO)
all_GO <- all_GO %>% replace(is.na(.), "n.s")
# save
write.csv(all_GO, "results/GO/All_GOenrichment.csv", row.names = F)

##### 4. Directional expression patterns for pre-exposure specific and infection specific genes
# get expression values, separate over/under expressed and pre-exposure specific and Bd infection specific
PdatO <- assay(vst$toe)[res_preExposure_specific[res_preExposure_specific$sigOvrexp_comb == "0011","gene"], ]
PdatU <- assay(vst$toe)[res_preExposure_specific[res_preExposure_specific$sigUndexp_comb == "0011","gene"], ]
BdatO <- assay(vst$toe)[res_BdExp2_specific[res_BdExp2_specific$sigOvrexp_comb == "1111","gene"], ]
BdatU <- assay(vst$toe)[res_BdExp2_specific[res_BdExp2_specific$sigUndexp_comb == "1111","gene"], ]
# scale expression per gene
PdatO <- t(scale(t(PdatO)))
PdatU <- t(scale(t(PdatU)))
BdatO <- t(scale(t(BdatO)))
BdatU <- t(scale(t(BdatU)))
# extract subsets
mean_scaled_exp <- list()
mean_scaled_exp$Pre_OE_TP1_mock.mock <- PdatO[,grep("UVC-TP1-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP1_bd.mock <- PdatO[,grep("\\dVC-TP1-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP1_mock.bd <- PdatO[,grep("UV-TP1-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP1_bd.bd <- PdatO[,grep("\\dV-TP1-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP2_mock.mock <- PdatO[,grep("UVC-TP2-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP2_bd.mock <- PdatO[,grep("\\dVC-TP2-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP2_mock.bd <- PdatO[,grep("UV-TP2-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_OE_TP2_bd.bd <- PdatO[,grep("\\dV-TP2-TOE", colnames(PdatO))]
mean_scaled_exp$Pre_UE_TP1_mock.mock <- PdatU[,grep("UVC-TP1-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP1_bd.mock <- PdatU[,grep("\\dVC-TP1-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP1_mock.bd <- PdatU[,grep("UV-TP1-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP1_bd.bd <- PdatU[,grep("\\dV-TP1-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP2_mock.mock <- PdatU[,grep("UVC-TP2-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP2_bd.mock <- PdatU[,grep("\\dVC-TP2-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP2_mock.bd <- PdatU[,grep("UV-TP2-TOE", colnames(PdatU))]
mean_scaled_exp$Pre_UE_TP2_bd.bd <- PdatU[,grep("\\dV-TP2-TOE", colnames(PdatU))]
mean_scaled_exp$BdI_OE_TP1_mock.mock <- BdatO[,grep("UVC-TP1-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP1_bd.mock <- BdatO[,grep("\\dVC-TP1-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP1_mock.bd <- BdatO[,grep("UV-TP1-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP1_bd.bd <- BdatO[,grep("\\dV-TP1-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP2_mock.mock <- BdatO[,grep("UVC-TP2-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP2_bd.mock <- BdatO[,grep("\\dVC-TP2-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP2_mock.bd <- BdatO[,grep("UV-TP2-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_OE_TP2_bd.bd <- BdatO[,grep("\\dV-TP2-TOE", colnames(BdatO))]
mean_scaled_exp$BdI_UE_TP1_mock.mock <- BdatU[,grep("UVC-TP1-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP1_bd.mock <- BdatU[,grep("\\dVC-TP1-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP1_mock.bd <- BdatU[,grep("UV-TP1-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP1_bd.bd <- BdatU[,grep("\\dV-TP1-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP2_mock.mock <- BdatU[,grep("UVC-TP2-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP2_bd.mock <- BdatU[,grep("\\dVC-TP2-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP2_mock.bd <- BdatU[,grep("UV-TP2-TOE", colnames(BdatU))]
mean_scaled_exp$BdI_UE_TP2_bd.bd <- BdatU[,grep("\\dV-TP2-TOE", colnames(BdatU))]
# get mean scaled expression for each gene in each group
for(i in names(mean_scaled_exp)){
  mean_scaled_exp[[i]] <- rowMeans(mean_scaled_exp[[i]])
}
# get data for ggplot
ggdat <- list()
# Pre_OE
len <- length(mean_scaled_exp$Pre_OE_TP1_mock.mock)
ggdat$Pre_OE_TP1 <- data.frame(expression = c(mean_scaled_exp$Pre_OE_TP1_mock.mock, mean_scaled_exp$Pre_OE_TP1_bd.mock, mean_scaled_exp$Pre_OE_TP1_mock.bd, mean_scaled_exp$Pre_OE_TP1_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
ggdat$Pre_OE_TP2 <- data.frame(expression = c(mean_scaled_exp$Pre_OE_TP2_mock.mock, mean_scaled_exp$Pre_OE_TP2_bd.mock, mean_scaled_exp$Pre_OE_TP2_mock.bd, mean_scaled_exp$Pre_OE_TP2_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
# Pre_UE
len <- length(mean_scaled_exp$Pre_UE_TP1_mock.mock)
ggdat$Pre_UE_TP1 <- data.frame(expression = c(mean_scaled_exp$Pre_UE_TP1_mock.mock, mean_scaled_exp$Pre_UE_TP1_bd.mock, mean_scaled_exp$Pre_UE_TP1_mock.bd, mean_scaled_exp$Pre_UE_TP1_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
ggdat$Pre_UE_TP2 <- data.frame(expression = c(mean_scaled_exp$Pre_UE_TP2_mock.mock, mean_scaled_exp$Pre_UE_TP2_bd.mock, mean_scaled_exp$Pre_UE_TP2_mock.bd, mean_scaled_exp$Pre_UE_TP2_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
# BdI_OE
len <- length(mean_scaled_exp$BdI_OE_TP1_mock.mock)
ggdat$BdI_OE_TP1 <- data.frame(expression = c(mean_scaled_exp$BdI_OE_TP1_mock.mock, mean_scaled_exp$BdI_OE_TP1_bd.mock, mean_scaled_exp$BdI_OE_TP1_mock.bd, mean_scaled_exp$BdI_OE_TP1_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
ggdat$BdI_OE_TP2 <- data.frame(expression = c(mean_scaled_exp$BdI_OE_TP2_mock.mock, mean_scaled_exp$BdI_OE_TP2_bd.mock, mean_scaled_exp$BdI_OE_TP2_mock.bd, mean_scaled_exp$BdI_OE_TP2_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
# BdI_UE
len <- length(mean_scaled_exp$BdI_UE_TP1_mock.mock)
ggdat$BdI_UE_TP1 <- data.frame(expression = c(mean_scaled_exp$BdI_UE_TP1_mock.mock, mean_scaled_exp$BdI_UE_TP1_bd.mock, mean_scaled_exp$BdI_UE_TP1_mock.bd, mean_scaled_exp$BdI_UE_TP1_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
ggdat$BdI_UE_TP2 <- data.frame(expression = c(mean_scaled_exp$BdI_UE_TP2_mock.mock, mean_scaled_exp$BdI_UE_TP2_bd.mock, mean_scaled_exp$BdI_UE_TP2_mock.bd, mean_scaled_exp$BdI_UE_TP2_bd.bd), group = c(rep("mock-mock", len), rep("Bd-mock", len), rep("mock-Bd", len), rep("Bd-Bd", len)))
rm(len)
# run paired wilcoxon tests for comparisons of interest
t_test <- list()
for(i in names(ggdat)){
  my.tt <- wilcox_test(data = ggdat[[i]], 
                  formula = expression ~ group, 
                  comparisons = list(c("mock-mock", "Bd-mock"), c("mock-Bd", "Bd-Bd"), c("mock-mock", "mock-Bd"), c("Bd-mock", "Bd-Bd")), 
                  p.adjust.method = "bonferroni", paired = TRUE
                  )
  my.tt <- my.tt %>% add_y_position()
  t_test[[i]] <- my.tt
  rm(my.tt)
}
# plot
# Pre-exposure specific 
# TP1 over-expression
p1 <- ggboxplot(ggdat$Pre_OE_TP1, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$Pre_OE_TP1, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(a) TP1 over-expressed genes") +
  theme(legend.position = "none")
# TP2 over-expression
p2 <- ggboxplot(ggdat$Pre_OE_TP2, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$Pre_OE_TP2, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(b) TP2 over-expressed genes") +
  theme(legend.position = "none")
# TP1 under-expression
p3 <- ggboxplot(ggdat$Pre_UE_TP1, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$Pre_UE_TP1, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(c) TP1 under-expressed genes") +
  theme(legend.position = "none")
# TP2 under-expression
p4 <- ggboxplot(ggdat$Pre_UE_TP2, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$Pre_UE_TP2, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(d) TP2 under-expressed genes") +
  theme(legend.position = "none")
plot_grid(p1, p2, p3, p4, ncol = 2)

# infection specific 
# TP1 over-expression
p1 <- ggboxplot(ggdat$BdI_OE_TP1, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$BdI_OE_TP1, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(a) TP1 over-expressed genes") +
  theme(legend.position = "none")
# TP2 over-expression
p2 <- ggboxplot(ggdat$BdI_OE_TP2, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$BdI_OE_TP2, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(b) TP2 over-expressed genes") +
  theme(legend.position = "none")
# TP1 under-expression
p3 <- ggboxplot(ggdat$BdI_UE_TP1, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$BdI_UE_TP1, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(c) TP1 under-expressed genes") +
  theme(legend.position = "none")
# TP2 under-expression
p4 <- ggboxplot(ggdat$BdI_UE_TP2, x = "group", y = "expression", color = "group", add = "jitter", shape = "group", add.params = list(size = 2, alpha = 0.5)) +
  scale_color_manual(values = c("mock-mock" = "#b38711", "mock-Bd" = "#b38711", "Bd-mock" = "#732f30", "Bd-Bd" = "#732f30")) +
  scale_shape_manual(values = c("mock-mock" = 19, "mock-Bd" = 17, "Bd-mock" = 19, "Bd-Bd" = 17)) +
  stat_pvalue_manual(t_test$BdI_UE_TP2, label = "p.adj.signif", hide.ns = F) +
  ylab("Mean scaled expression") + 
  xlab("Group") + 
  ggtitle("(d) TP2 under-expressed genes") +
  theme(legend.position = "none")
plot_grid(p1, p2, p3, p4, ncol = 2)