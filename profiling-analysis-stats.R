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
