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
# from = rep("utr-analysis", nrow(pre.gtf.df))
# score = rep(0.0, nrow(pre.gtf.df))
# phase = rep(".", nrow(pre.gtf.df))
write.table(pre.gtf.df, 
            file = "pre-gtf",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
#
# FINAL FORMATTING VIA AWK
# GETTING THE FORMATTING RIGHT
# awk 'BEGIN{ OFS = "\t" } { print $1, "utr-analysis", $6, $2, $3, "0.0", $5, ".", $7 " " $8 " " $9 " " $10}' pre-gtf | head
#
# THE ACTUAL PRINTING
# awk 'BEGIN{ OFS = "\t" } { print $1, "utr-analysis", $6, $2, $3, "0.0", $5, ".", $7 " " $8 " " $9 " " $10}' pre-gtf > sac_cer_yassour_utr.gtf
