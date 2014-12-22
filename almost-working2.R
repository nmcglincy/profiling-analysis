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
gtf.l.utrs = mclapply(gtf.l, UtrCreator, mc.cores = 4)
# 
# A GRAPH FOR A BIT OF A SANITY CHECK; MIGHT A TABLE BE BETTER
png("sanity-check.png")
plot(elementLengths(gtf.l), elementLengths(gtf.l.utrs))
abline(a = 0, b = 1, col = "red")
dev.off()
