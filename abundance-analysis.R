vanilla.data = read.delim("mTIF-anno-choco-ypd1",
                          header = FALSE,
                          sep = "\t",
                          stringsAsFactors = FALSE)
head(vanilla.data)
colnames(vanilla.data) = c("chr", "strand", "start", "end", "reads", "mTIF.class", "gene.id")
sum(is.na(vanilla.data$gene.id))
# [1] 6279
length(unique(vanilla.data$gene.id[!is.na(vanilla.data$gene.id)]))
# [1] 15149 - I think this must come from all the intergenic transcript names I made
# 
# I think I'm just going to take the covering-one-intact-orf class to make things a bit easier
vanilla.coio = subset(vanilla.data, mTIF.class == "Covering_one_intact_ORF")
head(vanilla.coio)
sum(is.na(vanilla.coio$gene.id))
# [1] 6279
length(unique(vanilla.coio$gene.id[!is.na(vanilla.coio$gene.id)]))
# [1] 4729
length(vanilla.coio$gene.id)
# [1] 142966
# 
library(ggplot2)
ggplot(vanilla.coio,
       aes(x = gene.id, y = log10(reads))) +
  geom_boxplot() +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14),
        axis.text.x = element_blank(),
        axis.text.y  = element_text(size=12),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

