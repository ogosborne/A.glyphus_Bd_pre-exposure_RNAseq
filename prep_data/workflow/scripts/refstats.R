# snakemake vars
files <- snakemake@input
out.tab <- snakemake@output[["tab"]]
out.fig <- snakemake@output[["fig"]]
# read files
inputs <- list()
for(n in seq_along(files)){
  # get comment lines
  file <- files[[n]]
  comm <- readLines(file)
  comm <- comm[substr(comm, 1, 1) == "#"]
  # get name
  nam <- gsub(".refstats", "", file)
  # get N unmapped
  reads <- as.numeric(gsub("[^0-9.-]", "", comm[grep("#Reads", comm)]))
  mapped <- as.numeric(gsub("[^0-9.-]", "", comm[grep("#Mapped", comm)]))
  unmapped <- data.frame(Ref = "unmapped", N.frag =(reads-mapped)/2)
  colnames(unmapped) <- c("Ref", nam)
  # get table
  tab <- read.table(file, comment.char = "#", header = F)[,c(1,8)]
  colnames(tab) <- c("Ref", nam) 
  tab <- rbind(unmapped, tab)
  inputs[[nam]] <- tab 
}
# combine tables
comb.tab <- Reduce(function(x, y) merge(x, y, by="Ref", sort = F), inputs)
comb.tab.m <- comb.tab
rownames(comb.tab.m) <- comb.tab.m$Ref
comb.tab.m <- as.matrix(comb.tab.m[,2:ncol(comb.tab.m)])
# barplot
pdf(out.fig)
layout(matrix(c(1,1,1,2),ncol = 4))
comb.tab.m.M <- comb.tab.m/1000000
barplot(comb.tab.m.M, ylab = "N Million Fragments", col = rainbow(nrow(comb.tab.m)), las = 3, cex.names = 0.4)
plot.new()
legend("center", legend = rownames(comb.tab.m),box.lty = 0, fill = rainbow(nrow(comb.tab.m)))
# histograms
layout(matrix(c(1,2,3,4),ncol = 2))
for(i in rownames(comb.tab.m)){
  hist(matrix(comb.tab.m[i,]), main = i, prob = TRUE, xlab = "N Fragments", col="peachpuff")
  lines(density(matrix(comb.tab.m[i,])), lwd = 2, col = "chocolate3")
}
dev.off()
# save table
write.csv(t(comb.tab.m), file = out.tab, row.names = T) # nolint: T_and_F_symbol_linter.
