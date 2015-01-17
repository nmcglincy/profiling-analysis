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
#   xlim(c(-100, 60)) +
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
# head(start.pos.df2)
# str(start.pos.df2)
# summary(start.pos.df2)
# 
# MAKING A LONG VERSION OF DF2
start.pos.ldf2 = start.pos.df2 %>%
  gather(reads, count, total.reads:reads.l34)
head(start.pos.ldf2)
ldf2.nototal = start.pos.ldf2 %>%
  filter(reads != "total.reads")
head(ldf2.nototal)
ggplot(ldf2.nototal, aes(x = pos, y = count, fill = reads)) +
  geom_bar(stat = "identity",
           position = "fill") +
  facet_wrap(~.id, ncol = 2) +
  xlab(expression(paste("Position relative to ", underline(bold(A)), "UG", sep = ""))) +
  ylab("Proportion of reads") +
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
# YEAH, DOESN'T LOOK THAT USEFUL
ggplot(ldf2.nototal, aes(x = pos, y = reads, fill = count)) +
  geom_tile() +
  scale_fill_continuous(na.value = "transparent", low = "white", high = "purple") +
  facet_wrap(~.id, ncol = 2) +
  xlab(expression(paste("Position relative to ", underline(bold(A)), "UG", sep = ""))) +
  ylab("Read length") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
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
# DOESN'T SEEM LIKE LIKE THERE'S MUCH DIFFERENCE
# APPEARANCE OF LONGER READS IN THE CDS - NI & LL SUGGES COMPARING RATIO OF NORMAL TO LONG READS
# BETWEEN UTR AND CDS TO SEE WHETHER IT'S A SAMPLING THING OR NOT
# MIGHT ALSO BE WORTH RELAXING LOWER LIMIT ON FP-FRAMING/COUNT TO SEE WHAT ELSE HAS COME THROUGH
# COULD ALSO JUST BE AN CONTAMINANT, LIKE THE +77 PEAK.
# 
# QUESTIONS FOR NICK
# ANY NORMALISATION FOR LIBRARY SIZE? - NOPE
# ANY NORMALISATION FOR UTR LENGTH? - NOPE
# WTF IS THE MASSIVE PEAK AT POS +77 - NI: SOME CONTAMINANT, e.g. AN INTRONIC snoRNA OR tRNA

# 
# END SITES
end.files = list.files(pattern = "*end_pos_len.txt")
end.pos.l = list()
for (i in end.files) {
  end.pos.l[[i]] = assign(i,
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
names(end.pos.l) = c("BY4742.A", "BY4742.B", "D19A", "D19B")
str(end.pos.l)
# library(plyr)
end.pos.df = ldply(end.pos.l)
head(end.pos.df)
# dim(end.pos.df)
# library(dplyr)
# library(tidyr)
end.pos.ldf = end.pos.df %>%
  gather(reads, count, total.reads:reads.l34)
head(end.pos.ldf)
dim(end.pos.ldf)
# 
# HOW MANY TOTAL READS IN EACH SAMPLE
end.pos.df %>%
  group_by(.id) %>%
  summarise(total = sum(total.reads))
# 
# library(ggplot2)
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(end.pos.df, aes(x = pos, y = total.reads, colour = .id)) +
  geom_line() +
  ylab("Total reads") +
  xlab(expression(paste("Position relative to ", "UG", underline(bold(A)), sep = ""))) +
  xlim(c(-40, 40)) +
#   coord_trans(ytrans = "log2") +
#   scale_y_continuous(breaks = c(seq(from = 20, to = 80, by = 20),
#                                 seq(from = 100, to = 800, by = 200),
#                                 seq(from = 1000, to = 5000, by = 1000))) +
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

end.ldf.nototal = end.pos.ldf %>%
  filter(reads != "total.reads")
head(end.ldf.nototal)
ggplot(end.ldf.nototal, aes(x = pos, y = count, fill = reads)) +
  geom_bar(stat = "identity",
           position = "fill") +
  facet_wrap(~.id, ncol = 2) +
  xlab(expression(paste("Position relative to ", underline(bold(A)), "UG", sep = ""))) +
  ylab("Proportion of reads") +
  scale_x_continuous(expand=c(0,0), limits = c(-40, 25)) + 
  scale_y_continuous(expand=c(0,0)) +
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
# YEAH, DOESN'T LOOK THAT USEFUL
ggplot(end.ldf.nototal, aes(x = pos, y = reads, fill = count)) +
  geom_tile() +
  scale_fill_continuous(na.value = "transparent", low = "white", high = "purple3") +
  facet_wrap(~.id, ncol = 2) +
  xlab(expression(paste("Position relative to ", "UG", underline(bold(A)), sep = ""))) +
  ylab("Read length") +
  scale_x_continuous(expand=c(0,0), limits = c(-40, 25)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
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
# DID ANYONE EVER SEE THIS DENSITY OF 30nt reads at the end?
