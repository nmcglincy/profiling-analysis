library(rtracklayer)
library(GenomicRanges)
# 
# READ IN GTF FILE
gtf = import("sac_cer_yassour.gtf", format = "GFF", asRangedData = FALSE)
# 
# ASSIGN TRANSCRIPT STRUCTURE AS NAME
names(gtf) = mcols(gtf)$type
# 
# SPLIT INTO GRangesList BY GENE NAME
gtf.l = split(gtf, gtf$group)
# 
# LAPPLY UTR IDENTIFYING FUNTION
source("utr-creator2.R")
gtf.l.utrs = mclapply(gtf.l, UtrCreator, mc.cores = 2)
# 
# A GRAPH FOR A BIT OF A SANITY CHECK; MIGHT A TABLE BE BETTER
png("sanity-check.png")
plot(elementLengths(gtf.l), elementLengths(gtf.l.utrs))
abline(a = 0, b = 1, col = "red")
dev.off()
# 
# LOOK-SEE
library(plyr)
for (i in 1:length(gtf.l.utrs)) {
  mcols(gtf.l.utrs[[i]])[,"type"] = names(gtf.l.utrs[[i]])
  mcols(gtf.l.utrs[[i]])[,"gene"] = names(gtf.l.utrs)[[i]]
}
head(gtf.l.utrs)
gtf.l.utrs.df = lapply(gtf.l.utrs, as.data.frame, row.names = NULL)
head(gtf.l.utrs.df)
pre.gtf.df = ldply(gtf.l.utrs.df)
pre.gtf.df = pre.gtf.df[,-1]
head(pre.gtf.df)


