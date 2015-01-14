# R ANALYSIS OF FP-FRAMING RESULTS FROM DELTA-TMA19 YEAST, USING sac_cer_yassour.gtf & sac_cer_yassour.bed
# 
# PERFORM LOCALLY IN /Users/nmcglincy/Documents/computing/github/nmcglincy/profiling-analysis
# HAVE TO DO THIS FROM THE TERMINAL TO ENTER MY PWD
# scp -r mcglincy@compute1.ingolia-lab.org:/mnt/ingolialab/mcglincy/NINM001/Sample_NINM01_index1/split3/old-gtf-alignments/Statistics .
#
# CONFIRM TRANSFER
setwd("./Statistics")
getwd()
list.files()
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
dev.off()
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
str(qexpr.data.l)
names(qexpr.data.l) = c("112A","112B","D19A","D19B")
qexpr.data.df = ldply(qexpr.data.l)
# head(qexpr.data.df)
names(qexpr.data.df)[1] = "lib"
# 
summary(qexpr.data.df)
# WTF, SOMETHING AT -3
qexpr.data.df[which(qexpr.data.df$lnt.quant.reg < 0),]
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
ls()
str(qexpr.data.l)
# 
# CHECKING THE LENGTH, ORDER AND COMPOSITION OF THE GENE.ID COLUMNS IS THE SAME
lapply(qexpr.data.l, nrow)
for (i in 1:length(qexpr.data.l)) {
  print(identical(qexpr.data.l[[1]]$gene.id, qexpr.data.l[[i]]$gene.id))
}


