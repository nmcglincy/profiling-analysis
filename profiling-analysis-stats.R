# R ANALYSIS OF FP-FRAMING RESULTS FROM DELTA-TMA19 YEAST, USING sac_cer_yassour.gtf & sac_cer_yassour.bed
# 
# PERFORM LOCALLY IN /Users/nmcglincy/Documents/computing/github/nmcglincy/profiling-analysis
# HAVE TO DO THIS FROM THE TERMINAL TO ENTER MY PWD
# scp -r mcglincy@compute1.ingolia-lab.org:/mnt/ingolialab/mcglincy/NINM001/Sample_NINM01_index1/split3/old-gtf-alignments/Statistics .
#
# CONFIRM TRANSFER
# setwd("./Statistics")
# getwd()
# list.files()
# 
# ANALYSIS OF FOOTPRINT LENGTH DISTRIBUTION
# READ IN FILES
frame.len.files = list.files(pattern = "*frame_len.txt")
frame.len.l = list()
for (i in frame.len.files) {
  frame.len.l[[i]] = assign(i, 
                            read.table(file = i, 
                                       header = FALSE,
                                       sep = "\t"))
}
names(frame.len.l) = c("112A", "112B", "D19A", "D19B")
frame.len.l = lapply(frame.len.l, function(x) {x = x[,1:2]})
# 
# REFORMAT INTO A DF
library(plyr)
frame.len.df = ldply(frame.len.l)
names(frame.len.df) = c("library", "rpf.length", "prop.reads")
# 
# PLOT
library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(frame.len.df, aes(x = rpf.length, y = prop.reads, colour = library)) +
  geom_line(size = 1.5) +
  ylab("Proportion of footprints") +
  xlab("Footprint length, nt") +
  ggtitle("Distribution of Footprint Lengths") +
  scale_x_continuous(breaks = seq(from = 25, to = 34, by = 1)) +
  scale_colour_manual(name = "Library", values = cbbPalette[c(1,2,3,8)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +    
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("dist-footprint-lengths.png", dpi = 400)

# 
# DIFFERENTIAL EXPRESSION ANALYSIS
# 
# READING IN THE DATA
qexpr.files = list.files(pattern = "*qexpr.txt")
qexpr.data.l = list()
for (i in qexpr.files) {
  qexpr.data.l[[i]] = assign(i,
                             read.table(file = i,
                              header = FALSE,
                              sep = "\t",
                              stringsAsFactors = FALSE))
  names(qexpr.data.l[[i]]) = c("gene.id", "lnt.quant.reg", "read.count")     
}
# str(qexpr.data.l)
names(qexpr.data.l) = c("112A","112B","D19A","D19B")
qexpr.data.df = ldply(qexpr.data.l)
# head(qexpr.data.df)
names(qexpr.data.df)[1] = "lib"
# 
# summary(qexpr.data.df)
# WTF, SOMETHING AT -3
# qexpr.data.df[which(qexpr.data.df$lnt.quant.reg < 0),]
# names(qexpr.data.df)
# 
# EDA
# RAW COUNTS
# HISTOGRAM BY LIB
ggplot(qexpr.data.df, aes(x = log2(read.count))) +
  geom_density(size = 1.25) +
  facet_wrap(~ lib, ncol = 2) +
  ylab("Density") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("read-count-density.png")

ggplot(qexpr.data.df, aes(x = log2(read.count), colour = lib)) +
  stat_ecdf() +
  ylab("ECDF") +
  scale_colour_manual(name = "Sample",
                      values = cbbPalette[c(1,2,3,8)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("read-count-ecdf.png")
# 
# THE D19 LIBS SEEM TO HAVE LESS READS PER GENE
# TOTAL READ.COUNT FOR EACH LIBRARY
library(dplyr)
write.csv(qexpr.data.df %>%
            group_by(lib) %>%
            summarise(total.reads = sum(read.count)),
          file = "total-reads-summ.csv",
          quote = FALSE,
          row.names = FALSE)
# 
# GRAPHS BASED ON READ.COUNT NORMALISED FOR TOTAL.READS
qexpr.data.df.nm = qexpr.data.df %>%
  group_by(lib) %>%
  mutate(norm.counts = read.count/sum(read.count))
# head(qexpr.data.df.nm)
# 
ggplot(qexpr.data.df.nm, aes(x = log2(norm.counts))) +
  geom_density(size = 1.25) +
  facet_wrap(~ lib, ncol = 2) +
  ylab("Density") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("norm-read-count-density.png")
ggplot(qexpr.data.df.nm, aes(x = log2(norm.counts), colour = lib)) +
  stat_ecdf() +
  ylab("ECDF") +
  scale_colour_manual(name = "Sample",
                      values = cbbPalette[c(1,2,3,8)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("norm-read-count-ecdf.png")
# THIS APPEARS TO REMOVE THE DIFFERENCE IN THE DISTRIBUTION BETWEEN THE 112S AND THE D19S
# SO PROBABLY THE PER GENE DECREASE IS JUST BECAUSE THERE ARE SLIGHTLY LESS READS IN THE D19 LIBS
# 
# DESEQ ANALYSIS
# CREATING A MATRIX OF COUNT DATA
# ls()
# str(qexpr.data.l)
# 
# CHECKING THE LENGTH, ORDER AND COMPOSITION OF THE GENE.ID COLUMNS IS THE SAME
# lapply(qexpr.data.l, nrow)
# for (i in 1:length(qexpr.data.l)) {
#   print(identical(qexpr.data.l[[1]]$gene.id, qexpr.data.l[[i]]$gene.id))
# }
# INDEED THEY ARE 
#
count.data.df = data.frame(row.names = qexpr.data.l[[1]]$gene.id,
  "112A" = qexpr.data.l[[1]]$read.count,
  "112B" = qexpr.data.l[[2]]$read.count,
  "D19A" = qexpr.data.l[[3]]$read.count,
  "D19B" = qexpr.data.l[[4]]$read.count)
# head(count.data.df)
count.data.m = as.matrix(count.data.df)
# class(count.data.m)
# head(count.data.m)
pheno.df = data.frame(row.names = colnames(count.data.m),
                      genotype = factor(c(rep("wt", 2), rep("del", 2)), levels = c("wt", "del")))
# str(pheno.df)
# colnames(count.data.m)
#
library(DESeq2)
deseq.data <- DESeqDataSetFromMatrix(countData = count.data.m,
                                     colData = pheno.df,
                                     design = ~ genotype)
deseq.data.ex <- DESeq(deseq.data)
res <- results(deseq.data.ex)
capture.output(summary(res), file = "result-summary.txt")
# SEEMS LIKE A LOT OF GENES EXCLUDED DUE TO LOW COUNTS
plotMA(res, ylim = c(-1.5,1.5))
# 
res.df = subset(as.data.frame(res), baseMean > 0 )
# 
# CREATE VARIABLE FOR SIGNIFICANT CHANGE
res.df = res.df %>%
  mutate(sig.change = padj < 0.1)
res.df$sig.change[is.na(res.df$sig.change)] = "indp.filt"
# 
# MY OWN MA PLOT
ggplot(res.df, aes(x = log10(baseMean), 
                   y = log2FoldChange,
                   colour = sig.change)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, size = 0.5) +
  xlab("log10(mean expression)") +
  scale_y_continuous(name = "log2(fold change)",
                     limits = c(-1.5, 1.5)) +
  scale_colour_manual(name = "padj < 0.1",
                      values = cbbPalette[c(6,1,2)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + 
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
  ggsave("ma-plot.png")
# 
# WRITE ORDERED RESULTS AS A CSV
resOrdered <- res[order(res$padj),]
write.csv(resOrdered,
          file = "res-ordered.csv",
          quote = FALSE)
sig.changes = subset(resOrdered, padj < 0.1)
sig.changes.ord = sig.changes[order(sig.changes$log2FoldChange),]
# 
# SOME SINGLE GENE PLOT: TMA19, PLUS THE BIGGEST CHANGE EITHER WAY
YKL056C.dat = plotCounts(deseq.data.ex, gene="YKL056C", intgroup="genotype", returnData = TRUE)
YBR093C.dat = plotCounts(deseq.data.ex, gene="YBR093C", intgroup="genotype", returnData = TRUE)
YHR216W.dat = plotCounts(deseq.data.ex, gene="YHR216W", intgroup="genotype", returnData = TRUE)

ggplot(YKL056C.dat,
       aes(x = genotype, y = count)) +
  geom_point(size = 5, position = position_jitter(w = 0.1, h = 0)) +
  coord_trans(ytrans = "log10") +
  scale_y_continuous(limits = c(1,4000),
    breaks = c(1, 10, 100, 1000, 10000), 
                     labels = c(1, 10, 100, 1000, 10000)) +
  ylab("Moderated read count") +
  xlab("Genotype") +
  scale_x_discrete(labels = c(expression(italic("BY4742")), 
                              expression(paste(Delta, italic("tma19"), sep = "")))) +
  ggtitle(expression(paste(italic("tma19/YKL056C"), " expression", sep = ""))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18, vjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("YKL056C-counts.png")

ggplot(YBR093C.dat,
       aes(x = genotype, y = count)) +
  geom_point(size = 5, position = position_jitter(w = 0.1, h = 0)) +
  coord_trans(ytrans = "log10") +
#   scale_y_continuous(limits = c(1,4000),
#                      breaks = c(1, 10, 100, 1000, 10000), 
#                      labels = c(1, 10, 100, 1000, 10000)) +
  ylab("Moderated read count") +
  xlab("Genotype") +
  scale_x_discrete(labels = c(expression(italic("BY4742")), 
                              expression(paste(Delta, italic("tma19"), sep = "")))) +
  ggtitle(expression(paste(italic("YBR093C"), " expression", sep = ""))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18, vjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("YBR093C-counts.png")

ggplot(YHR216W.dat,
       aes(x = genotype, y = count)) +
  geom_point(size = 5, position = position_jitter(w = 0.1, h = 0)) +
  coord_trans(ytrans = "log10") +
  #   scale_y_continuous(limits = c(1,4000),
  #                      breaks = c(1, 10, 100, 1000, 10000), 
  #                      labels = c(1, 10, 100, 1000, 10000)) +
  ylab("Moderated read count") +
  xlab("Genotype") +
  scale_x_discrete(labels = c(expression(italic("BY4742")), 
                              expression(paste(Delta, italic("tma19"), sep = "")))) +
  ggtitle(expression(paste(italic("YHR216W"), " expression", sep = ""))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 1, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18, vjust = 1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
ggsave("YHR216W-counts.png")
# 
# TRYING WITHOUT THE INDEPENDENT FILTERING
# res2 = results(deseq.data.ex,
#                independentFiltering = FALSE)
# summary(res2)
# res2Ordered = res2[order(res2$padj),]
# head(res2Ordered)
# plotMA(res2, ylim = c(-1.5,1.5))
# TAKING OUT THE INDEPENDANT FILTERED RESULTED IN LESS CHANGES CALLED AS SIGNIFICANT
# 
# DISPERSION ESTIMATES
png("dispersion-plot.png")
plotDispEsts(deseq.data.ex)
dev.off()
# 
# TRANSFORMED COUNTS HCLUST TOGETHER AS EXPECTED, USING CODE FROM VIGNETTE
rld <- rlog(deseq.data.ex)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
mat
rownames(mat) <- colnames(mat) <- with(colData(deseq.data.ex),
                                       paste(genotype))
hc <- hclust(distsRL)
# str(hc)
library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

png("rlog-heatmap.png")
heatmap.2(mat, 
          Rowv=as.dendrogram(hc),
          symm=TRUE, 
          trace="none",
          col = rev(hmcol), 
          margin=c(13, 13))
dev.off()
# 
# GO ANALYSIS
# 
# RATHER THAT JUST TAKE THE WHOLE YEAST GENOME AS MY BACKGROUND, I WANTED TO TAKE THE GENES THAT I
# HAD BEEN ABLE TO DETECT AT ALL
library(GO.db)
ls("package:GO.db")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Sc.sgd.db")
library("org.Sc.sgd.db")
ls("package:org.Sc.sgd.db")




