#
# SEPARATING OUT THE UTRS IN NICK'S YEAST GTF/BED
#
# TOPHAT ETC WANT A GTF, SO WE'LL TRY TO START FROM THAT.
# 
# OPTIONS FOR IMPORT
# rtracklayer:import(format = "GFF")
# GenomicRanges:makeGRangesFromDataFrame()
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
nick.gtf = import(con = "sac_cer_yassour.gtf",
                  format = "GFF")
head(nick.gtf)
names(nick.gtf) = nick.gtf$type
head(nick.gtf)
mcols(nick.gtf) = mcols(nick.gtf)$group
# 
#
nick.grl = split(nick.gtf, mcols(nick.gtf)$value)
#
# SOME SANITY CHECKS
head(nick.grl)
lapply(nick.grl, mcols)
hist(elementLengths(nick.grl))
length(nick.grl)
length(nick.gtf)
length(unique(mcols(nick.gtf)$group))
# 
# SPLITTING INTO UTRS
mcols(nick.grl[[1]]) = NULL

source("utr-creator.R")
nick.grl[1]
UtrCreator(nick.grl[1])
foo = lapply(nick.grl, UtrCreator)
foo = nick.grl[[1]]
foo
UtrCreator(foo)
disjoin(foo)
foo
c(foo, disjoin(foo))

mcols(foo) = NULL
foo
UtrCreator(foo)
# THE META DATA IS GETTING IN THE WAY
nick.grl1 = mcols(nick.grl) = NULL
head(nick.grl1)

mcols(nick.grl)

head(nick.grl)
bar = lapply(nick.grl, function(x) {elementMetadata(x) = NULL})
head(bar)

str(nick.grl)

bar = unlist(nick.grl)
head(bar)
mcols(bar) = NULL


bar = ranges(nick.grl, use.mcols = FALSE)
head(bar)
rm(list = ls())

length(nick.grl)
foo = endoapply(nick.grl, function(x) {elementMetadata(x)<-NULL; x})
bar = endoapply(nick.grl, function(x) {mcols(x)<-NULL; x})
head(foo)
head(bar)
?ranges

nick.grl1 = endoapply(bar, UtrCreator)
UtrCreator(bar[1])
