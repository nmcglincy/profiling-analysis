# Creating a meta-report of the fastx-split, rRNA alignments and genome alignments
# 
# To be run from within the split directory
# Result of fastx-split
split.result = read.table("fates.txt",
                          sep = "\t",
                          header = FALSE,
                          stringsAsFactors = FALSE,
                          col.names = c("sample",
                                        "index",
                                        "no.reads",
                                        "prop.library"))
# 
# rRNA alignment results
rrna.reports.files = list.files(pattern = "*rrna_count.txt")
rrna.reports.l = list()
for (i in rrna.reports.files) {
  rrna.reports.l[[i]] = assign(i, read.table(i,
                                             header = FALSE,
                                             sep = "\t",
                                             stringsAsFactors = FALSE,
                                             col.names = c("source",
                                                           "no.reads",
                                                           "prop.input"))) 
}
library(stringr)
names(rrna.reports.l) = str_sub(rrna.reports.files, 
                                start = 1,
                                end = str_locate(rrna.reports.files, pattern = "_")[,2]-1)
library(plyr)
library(dplyr)
rrna.reports.df = ldply(rrna.reports.l)
# 
# tophat alignment results
# prep_reads.info
prep.reads.files = str_c(names(rrna.reports.l), "_norrna_v_genome/prep_reads.info", sep = "")
prep.reads.l = list()
for (i in prep.reads.files) {
  prep.reads.l[[i]] = assign(i, read.table(i,
                                           header = FALSE,
                                           skip = 2,
                                           sep = "=",
                                           stringsAsFactors = FALSE,
                                           col.names = c("class", "no.reads")))
}
names(prep.reads.l) = names(rrna.reports.l)
prep.reads.df = ldply(prep.reads.l)
# 
# align_summary.txt
aln.summ.files = str_c(names(rrna.reports.l), "_norrna_v_genome/align_summary.txt", sep = "")
aln.summ.names = names(rrna.reports.l)
# 
aln.summ.l = list()
for (i in aln.summ.files) {
  aln.summ.l[[i]] = assign(i,
                           scan(file = i,
                                what = "character",
                                skip = 1,
                                nlines = 3,
                                strip.white = TRUE,
                                quiet = TRUE,
                                sep = c("\t", ":", "(", ")", "%", " ")))
  aln.summ.l[[i]] = as.numeric(c(unlist(str_extract_all(aln.summ.l[[i]], '[0-9]+$')), 
                                 unlist(str_extract_all(aln.summ.l[[i]], '[0-9]+ '))))
  aln.summ.l[[i]] = data.frame(var = c("tophat.input", "tophat.mapped", "tophat.om.multimapped", "tophat.omm.gt20"),
                               no.reads = aln.summ.l[[i]])
}
names(aln.summ.l) = aln.summ.names
aln.summ.df = ldply(aln.summ.l)
# 
# Putting all the files together
# 
# aln.summ.df
# prep.reads.df
# rrna.reports.df
# split.result
library(tidyr)
prep.reads.df2 = prep.reads.df %>% spread(class, no.reads)
rrna.reports.df2 = rrna.reports.df[,1:3] %>% spread(source, no.reads)
aln.summ.df2 = aln.summ.df %>% spread(var, no.reads)
foo = merge(split.result, 
            rrna.reports.df2,
            by.x = "sample", 
            by.y = ".id", 
            all = TRUE)
bar = merge(foo,
            prep.reads.df2,
            by.x = "sample",
            by.y = ".id",
            all = TRUE)
monkey = merge(bar,
               aln.summ.df2,
               by.x = "sample",
               by.y = ".id",
               all = TRUE)
banana = as.data.frame(t(monkey))
colnames(banana) = monkey$sample
ape = banana[-1,]
row.names(ape) = c("index",
                   "no.reads",
                   "prop.lib",
                   "rRNA.dep.input",
                   "non.rRNA.reads",
                   "rRNA.read",
                   "prep_info.reads.in",
                   "prep_info.reads.out",
                   "tophat.input",
                   "tophat.mapped",
                   "tophat.om.multimapped",
                   "tophat.omm.gt20")
ape = ape[c(1:4,6,5,7:12),]
ape = data.frame(var = row.names(ape),
                 ape)
# write ape out
write.table(ape,
            file = "split-aln-meta-report.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
# 