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
