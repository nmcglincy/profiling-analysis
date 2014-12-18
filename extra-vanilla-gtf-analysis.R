library(GenomicFeatures)
exVanilla = makeTranscriptDbFromGFF(file = "mTIF-exVanilla.gtf", 
                                    format = "gtf")
exVanilla
head(exVanilla)
transcripts(exVanilla)
metadata(exVanilla)
seqinfo(exVanilla)
str(as.list(exVanilla))
