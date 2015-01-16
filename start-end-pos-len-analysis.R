# ANALYSIS OF THE _start_pos_len.txt & _end_pos_len.txt FILES FROM FP-COUNT OUTPUT
# 20150116
# 
# START SITES
start.files = list.files(pattern = "*start_pos_len.txt")
start.pos.l = list()
for (i in start.files) {
  start.pos.l[[i]] = assign(i,
                            read.table(i,
                                       header = FALSE,
                                       skip = 2,
                                       sep = "\t",
                                       col.names = c("pos",
                                                     "total.reads",
                                                     paste("reads.l",
                                                           25:34,
                                                           sep = "")),
                                       stringsAsFactors = FALSE))
}
names(start.pos.l) = c("BY4742.A", "BY4742.B", "D19A", "D19B")
# str(start.pos.l)
library(plyr)
start.pos.df = ldply(start.pos.l)
# head(start.pos.df)
# dim(start.pos.df)
library(dplyr)
library(tidyr)
start.pos.ldf = start.pos.df %>%
  gather(reads, count, total.reads:reads.l34)
# head(start.pos.ldf)
# dim(start.pos.ldf)
# 
# HOW MANY TOTAL READS IN EACH SAMPLE
start.pos.df %>%
  group_by(.id) %>%
  summarise(total = sum(total.reads))
# 
library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(start.pos.df, aes(x = pos, y = total.reads, colour = .id)) +
  geom_line() +
  ylab("Total reads") +
  xlab(expression(paste("Position relative to ", underline(bold(A)), "UG", sep = ""))) +
  xlim(c(-100, 60)) +
  coord_trans(limy = c(2, 5000), ytrans = "log2") +
  scale_y_continuous(breaks = c(seq(from = 20, to = 80, by = 20),
                                seq(from = 100, to = 800, by = 200),
                                seq(from = 1000, to = 5000, by = 1000))) +
  scale_colour_manual(name = "Library", values = cbbPalette[c(1,2,3,8)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +    
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
# 
# SUBSET TO AVOID +77 PEAK PERMENANTLY
start.pos.df2 = subset(start.pos.df, pos < 77)
# start.pos.df2
head(start.pos.df2)
str(start.pos.df2)
summary(start.pos.df2)
ggplot(start.pos.df2, aes(x = pos, y = total.reads, colour = .id)) +
  stat_ecdf() +
  ylab("ECDF") +
  xlab(expression(paste("Position relative to ", underline(bold(A)), "UG", sep = ""))) +
  scale_colour_manual(name = "Library", values = cbbPalette[c(1,2,3,8)]) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +    
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y  = element_text(size=14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))


# QUESTIONS FOR NICK
# ANY NORMALISATION FOR LIBRARY SIZE?
# ANY NORMALISATION FOR UTR LENGTH?
# WTF IS THE MASSIVE PEAK AT POS +77 - NI: SOME CONTAMINANT, e.g. AN INTRONIC snoRNA OR tRNA

# 
# END SITES
