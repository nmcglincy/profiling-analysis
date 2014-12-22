library(rtracklayer)
library(GenomicRanges)
gtf = import("sac_cer_yassour.gtf", format = "GFF", asRangedData = FALSE)
bed = import("sac_cer_yassour.bed", format = "BED", asRangedData = FALSE)
head(gtf)
head(bed)
gtf$type
names(gtf) 
names(gtf) = mcols(gtf)$type
gtf.l = split(gtf, gtf$group)
names(gtf.l)
head(gtf.l)

names(gtf.l[1])
names(gtf.l[[1]])
gtf.l[1]$type

append(gtf.l[[1]], disjoin(gtf.l[[1]]))

unique(c(gtf.l[[6690]], disjoin(gtf.l[[6690]]), ignore.mcols = TRUE))

list.files()
source("utr-creator.R")
fooTastic = lapply(gtf.l, UtrCreator)
#Error in `unsafe.names<-`(`*tmp*`, value = c("exon", "CDS", "3pUTR")) : 
#  too many names
foo = endoapply(gtf.l, UtrCreator)
UtrCreator(gtf.l[[1]])
gtf.l[["gene_id \"YPL038W-A\"; transcript_id \"YPL038W-A\";"]]
length(gtf.l)
UtrCreator(gtf.l[[6690]])
disjoin(gtf.l[[6690]])
bar = 

foo = endoapply

?GRangesList

foo = disjoin(gtf.l)
head(foo)

unlist(lapply(foo, elementLengths))




library(GenomicFeatures)
?makeTranscriptDbFromGFF