library(GenomicFeatures)
exVanilla = makeTranscriptDbFromGFF(file = "mTIF-exVanilla.gtf", 
                                    format = "gtf")
exVanilla
head(exVanilla)
transcripts(exVanilla)
metadata(exVanilla)
seqinfo(exVanilla)
str(as.list(exVanilla))
#
?import
library(rtracklayer)
xv.data = import.bed(con = "mTIF-exVanilla.bed",
                     asRangedData = FALSE,
                     genome = "saccCer3")
xv.data
?split
foo = split(xv.data, mcols(xv.data)$name)
foo
length(mcols(xv.data)$name)
metadata(xv.data)
head(foo)
foo["YAL002W"]
str(foo)
length(unique(names(foo)))
class(foo)
length(foo)
seqlengths(foo)
bar = elementLengths(foo)
hist(bar)
