library(GenomicRanges)
library(GenomicFeatures)
ir1 = IRanges(start = c(1, 4), end = c(10, 8))
ir1
?disjoin
names(ir1) = c("exon", "CDS")
disjoin(ir1)
c(ir1, disjoin(ir1))
?unique
ir2 = unique(c(ir1, disjoin(ir1)))
ir2
names(ir2) = c(names(ir1), "5pUTR", "3pUTR")

UtrCreator = function(x) {
	require(GenomicRanges)
	require(GenomicFeatures)
	if (length(disjoin(x)) > length(x)) {
		x = unique(c(x, disjoin(x)))
	} else {
		x = x
	}
	names(x) = c("exon", "CDS", "5pUTR", "3pUTR")[1:length(x)]
	return(x)
}

goi1 = IRanges(start = c(1, 4), 
			   end = c(10, 8),
			   names = c("exon", "CDS"))
goi2 = IRanges(start = c(1, 2), 
			   end = c(15, 10),
			   names = c("exon", "CDS"))
gene.list = IRangesList(goi1 = goi1, goi2 = goi2)
gene.list

gene.list2 = lapply(gene.list, UtrCreator)
gene.list2
UtrCreator(goi1)
UtrCreator(goi2)

# NOW LET'S TRY IT FOR GRanges & GRangesList FORMATS

?GRanges
goi1
gr.go1 = GRanges(start = c(1, 4), 
			   end = c(10, 8),
			   names = c("exon", "CDS"), 
				 seqnames = c("chr1", "chr1"),
				 stand = Rle(c("+", "+")))
gr.go1				 

gr1 = GRanges(start = c(1, 4),
				 end = c(10, 8))
gr1 = Granges(seqname)

# I DON'T UNDERSTAND RLE				 
x = rev(rep(6:10, 1:5))
x
rle(x)
y = "WWWWWWWWWWWWBWWWWWWWWWWWWBBBWWWWWWWWWWWWWWWWWWWWWWWWBWWWWWWWWWWWWWW"
library(stringr)
y = unlist(str_split(y ,""))
y = y[[1]]
y
rle(y)
#
# NOW I UNDERSTAND RLE

> gr <-
	GRanges(seqnames =
	Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
	ranges =
 	IRanges(1:10, end = 7:16, names = head(letters, 10)),
 	strand =
 	Rle(strand(c("-", "+", "*", "+", "-")),
 	c(1, 2, 2, 3, 2)),
 	score = 1:10,
 	GC = seq(1, 0, length=10))

gr.co1 = GRanges(seqnames = Rle(c("chr1"), c(2)),
				 ranges = IRanges(start = c(1, 1), 
			   					  end = c(10, 8),
			   					  names = c("exon", "CDS")),
			   	 strand = Rle(strand("+")))
gr.co1
UtrCreator(gr.co1)
# THROWS AN ERROR BECAUSE YOU'RE TRYING TO FOR 4 NAMES ON 2 THINGS
length(gr.co1)
# COULD YOU USE THE LENGTH TO SAY HOW MANY NAMES
UtrCreator(gr.co1)
# NOW IT REMOVES THE CDS ANNOTATION, NEED IT TO TEST THE LENGTH OF THE DISJOIN
length(disjoin(gr.co1))
length(gr.co1)

if (length(disjoin(x)) > length(x)) {
	x = unique(c(x, disjoin(x)))
} else {
	x = x
}
# OK, NOW TEST IT
UtrCreator(gr.co1)
UtrCreator(gr.co2)
lapply(gr.list, UtrCreator)

gr.co2 = GRanges(seqnames = Rle(c("chr2"), c(2)),
				 ranges = IRanges(start = c(1, 2), 
			   					  end = c(15, 10),
				    			  names = c("exon", "CDS")),
			   	 strand = Rle(strand("-")))
UtrCreator(gr.co2)

gr.list = GRangesList(gr.co1, gr.co2)
gr.list
lapply(gr.list, UtrCreator)

#
# SCROLLING THROUGH NICK'S GTF HAS ILLUSTATED SOME EDGE CASES, E.G.
>1
chr03	bed-to-gtf	exon	107023	107033	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
chr03	bed-to-gtf	exon	107111	107191	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
chr03	bed-to-gtf	exon	107288	107521	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
chr03	bed-to-gtf	CDS	107023	107033	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
chr03	bed-to-gtf	CDS	107111	107191	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
chr03	bed-to-gtf	CDS	107288	107417	0.0	+	.	gene_id "YCL005W-A"; transcript_id "YCL005W-A"; 
>2
chr06	bed-to-gtf	exon	63974	64222	0.0	-	.	gene_id "YFL034C-B"; transcript_id "YFL034C-B"; 
chr06	bed-to-gtf	exon	62862	63859	0.0	-	.	gene_id "YFL034C-B"; transcript_id "YFL034C-B"; 
chr06	bed-to-gtf	CDS	63974	63993	0.0	-	.	gene_id "YFL034C-B"; transcript_id "YFL034C-B"; 
chr06	bed-to-gtf	CDS	63016	63859	0.0	-	.	gene_id "YFL034C-B"; transcript_id "YFL034C-B"; 

# DOESN'T WORK ON ASYMMETRIC CASES, NEED TO FIND SOME FUNCTION TO ASSESS THIS
foo = unique(c(gr.co1, disjoin(gr.co1)))

gr.co1 = GRanges(seqnames = Rle(c("chr1"), c(2)),
                 ranges = IRanges(start = c(1, 1), 
                                  end = c(10, 8),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("+")))
gr.co2 = GRanges(seqnames = Rle(c("chr2"), c(2)),
                 ranges = IRanges(start = c(1, 3), 
                                  end = c(10, 8),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("+")))
gr.co3 = GRanges(seqnames = Rle(c("chr3"), c(2)),
                 ranges = IRanges(start = c(1, 1), 
                                  end = c(10, 10),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("+")))
gr.co4 = GRanges(seqnames = Rle(c("chr1"), c(2)),
                 ranges = IRanges(start = c(1, 3), 
                                  end = c(10, 10),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("+")))



gr.list = GRangesList(gr.co1, gr.co2, gr.co3, gr.co4)
gr.list

source("utr-creator.R")

UtrCreator(gr.co1)
UtrCreator(gr.co2)
UtrCreator(gr.co3)
UtrCreator(gr.co4)
lapply(gr.list, UtrCreator)

start(gr.co1)
disjoin(gr.co1)
start(disjoin(gr.co1))
strand(gr.co1) == "+"
?nearest
disjoin(gr.co1)[2]
precede(disjoin(gr.co1)[2], disjoin(gr.co1)[1]) == 1
follow(disjoin(gr.co1)[2], disjoin(gr.co1)[1]) == NA

if (follow(disjoin(gr.co4)[2], disjoin(gr.co4)[1]) == 1)
gr.co4["CDS"]
setdiff(disjoin(gr.co4), gr.co4["CDS"])
setdiff(disjoin(gr.co1), gr.co1["CDS"])
precede( setdiff(disjoin(gr.co4), gr.co4["CDS"]), gr.co4["CDS"] ) == 1
follow( setdiff(disjoin(gr.co4), gr.co4["CDS"]), gr.co4["CDS"] ) == 1
foo = precede( setdiff(disjoin(gr.co1), gr.co1["CDS"]), gr.co1["CDS"] )
foo = as.character(foo)
class(foo)

as.character(precede( setdiff(disjoin(gr.co1), gr.co1["CDS"]), gr.co1["CDS"] )) == "NA"
foo = as.character(precede( setdiff(disjoin(gr.co1), gr.co1["CDS"]), gr.co1["CDS"] ))
is.na(foo)
foo = as.character(precede( setdiff(disjoin(gr.co1), gr.co1["CDS"]), gr.co1["CDS"] ))
is.na(precede( setdiff(disjoin(gr.co1), gr.co1["CDS"]), gr.co1["CDS"] ))
# [1] TRUE
is.na(precede( setdiff(disjoin(gr.co4), gr.co4["CDS"]), gr.co4["CDS"] ))
# [1] FALSE

# 
# TWO MORE EDGE CASES I CAN THINK OF:
# 1. HOW WILL PRECEDE/FOLLOW WORK ON THE NEGATIVE STRAND
# 2. WILL IT BE CONFUSED BY MULTI-EXON GENES
# 
# 1
gr.co5 = GRanges(seqnames = Rle(c("chr1"), c(2)),
                 ranges = IRanges(start = c(1, 1), 
                                  end = c(10, 8),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("-")))
gr.co6 = GRanges(seqnames = Rle(c("chr2"), c(2)),
                 ranges = IRanges(start = c(1, 3), 
                                  end = c(10, 8),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("-")))
gr.co7 = GRanges(seqnames = Rle(c("chr3"), c(2)),
                 ranges = IRanges(start = c(1, 1), 
                                  end = c(10, 10),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("-")))
gr.co8 = GRanges(seqnames = Rle(c("chr1"), c(2)),
                 ranges = IRanges(start = c(1, 3), 
                                  end = c(10, 10),
                                  names = c("exon", "CDS")),
                 strand = Rle(strand("-")))



gr.list = GRangesList(gr.co1, gr.co2, gr.co3, gr.co4, gr.co5, gr.co6, gr.co7, gr.co8, gr.co9)
gr.list

source("utr-creator.R")

UtrCreator(gr.co1)
UtrCreator(gr.co2)
UtrCreator(gr.co3)
UtrCreator(gr.co4)
lapply(gr.list, UtrCreator)
# 
# LOL, ACTUALLY THE NORMAL CASE WAS THE PROBLEM
strand(gr.co1) == "-"
unique(strand(gr.co1)) == "-"
# THAT FIX WORKS 
# 
# MULTI-EXON CASE
gr.co9 = GRanges(seqnames = Rle(c("chr2"), c(6)),
                 ranges = IRanges(start = c(10, 40, 70, 10, 40, 70), 
                                  end = c(30, 60, 100, 30, 60, 90),
                                  names = rep(c("exon", "CDS"), each = 3)),
                 strand = Rle(strand("+")))
gr.co9
UtrCreator(gr.co9)

disjoin(gr.co9)
gr.co9[names(gr.co9) == "CDS"]
gr.co9["CDS"]
names(gr.co9) == "CDS"
c(gr.co9, setdiff(disjoin(gr.co9), gr.co9[names(gr.co9) == "CDS"]))
