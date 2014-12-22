library(rtracklayer)
library(GenomicRanges)
gtf = import("sac_cer_yassour.gtf", format = "GFF", asRangedData = FALSE)
# bed = import("sac_cer_yassour.bed", format = "BED", asRangedData = FALSE)
head(gtf)
# head(bed)
# gtf$type
# names(gtf) 
names(gtf) = mcols(gtf)$type
gtf.l = split(gtf, gtf$group)
# names(gtf.l)
head(gtf.l)
# 
# names(gtf.l[1])
# names(gtf.l[[1]])
# gtf.l[1]$type
# 
# append(gtf.l[[1]], disjoin(gtf.l[[1]]))
# 
# unique(c(gtf.l[[6690]], disjoin(gtf.l[[6690]]), ignore.mcols = TRUE))
# 
# list.files()
source("utr-creator2.R")
fooTastic = lapply(gtf.l, UtrCreator)
#Error in `unsafe.names<-`(`*tmp*`, value = c("exon", "CDS", "3pUTR")) : 
#  too many names
length(gtf.l)
gtf.l[[6696]]
# THIS IS THE RIGHT SORT I THINK, HOW DOES IT WORK
UtrCreator(gtf.l[[6696]])
# THIS SEEMS TO WORK FINE, SO IT MUST BE FROM ONE OF THE MULTI-EXON CASES
which(elementLengths(disjoin(gtf.l)) > 3)
gtf.l[[6663]]
# THIS CASE HAS BOTH UTRS
UtrCreator(gtf.l[[6663]])
# THAT WORKED FINE
UtrCreator(gtf.l[[103]])
# THIS FAILS, BUT FOR A DIFFERENT REASON
disjoin(gtf.l[[103]])
# FIXED THIS, BUT IT'S STILL NOT WORKING
odd.ones = which(elementLengths(disjoin(gtf.l)) > 3)
is.vector(odd.ones)
for (i in 1:length(gtf.l)) {
  foo = UtrCreator(gtf.l[[i]])
  print(i)
  print(foo)
}
UtrCreator(gtf.l[[59]])
UtrCreator(gtf.l[[60]])
# NARROWED THE PROBLEM TO INSTANCE 60
gtf.l[[60]]

disjoin(gtf.l[[60]])
length(gtf.l[[60]])
length(disjoin(gtf.l[[60]])) == length(gtf.l[[60]])/2
# OK, THINK I FIXED THIS
UtrCreator(gtf.l[[60]])
# WORKS!
fooTastic = lapply(gtf.l, UtrCreator)
# SANITY CHECK GRAPH
plot(elementLengths(gtf.l), elementLengths(fooTastic))

UtrCreator(gtf.l[[5188]])
disjoin(gtf.l[[5188]])

endoapply(gtf.l, UtrCreator)

# UtrCreator(gtf.l[[1]])
# gtf.l[["gene_id \"YPL038W-A\"; transcript_id \"YPL038W-A\";"]]
# length(gtf.l)
# UtrCreator(gtf.l[[6690]])
# disjoin(gtf.l[[6690]])
# bar = 
# 
# foo = endoapply
# 
# ?GRangesList
# 
# foo = disjoin(gtf.l)
# head(foo)
# 
# unlist(lapply(foo, elementLengths))
# 
# 
# 
# 
# library(GenomicFeatures)
# ?makeTranscriptDbFromGFF
