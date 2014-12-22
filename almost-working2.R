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
head(gtf.l.utrs)
head(unlist(gtf.l.utrs))
foo = unlist(gtf.l.utrs)
head(foo)
class(foo)
class(gtf.l.utrs)
mcols(gtf.l.utrs)$gene.id = names(gtf.l.utrs)
names(gtf.l.utrs)
foo = lapply(gtf.l.utrs, as.data.frame)
names(gtf.l.utrs[[1]])
library(plyr)
foo = ldply(gtf.l.utrs)

foo = gtf.l.utrs[[1]]
foo
class(foo)
as.data.frame(foo)
seqnames(foo)
row.names(foo)
names(foo)
foo$type = names(foo)
foo
as.data.frame(foo, row.names = NULL)
bar = endoapply(gtf.l.utrs, function(x) {mcols(x)[,"type"] = names(x)})
head(bar)
str(bar)

length(foo)
length(gtf.l.utrs)

for (i in 1:length(gtf.l.utrs)) {
  mcols(gtf.l.utrs[[i]])[,"type"] = names(gtf.l.utrs[[i]])
  mcols(gtf.l.utrs[[i]])[,"gene"] = names(gtf.l.utrs)[[i]]
}

names(gtf.l.utrs)

mcols(gtf.l.utrs[[1]])[,"type"] = names(gtf.l.utrs[[1]])
mcols(gtf.l.utrs[[1]])[,"gene"] = names(gtf.l.utrs)[[1]]
# OK THESE WORK
gtf.l.utrs[[1]]
head(gtf.l.utrs)
