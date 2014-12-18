vanilla.data = read.delim("mTIF-anno-choco-ypd1",
                          header = FALSE,
                          sep = "\t",
                          stringsAsFactors = FALSE)
head(vanilla.data)
colnames(vanilla.data) = c("chr", "strand", "start", "end", "reads", "mTIF.class", "gene.id")
sum(is.na(vanilla.data$gene.id))
# [1] 6279
sum(unique(length(vanilla.data$gene.id)))
# [1] 238747 - I think this must come from all the intergenic transcript names I made
sum(unique(length(vanilla.data$gene.id))) - sum(is.na(vanilla.data$gene.id))
# [1] 232468
