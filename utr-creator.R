UtrCreator = function(x) {
  require(GenomicRanges)
  require(GenomicFeatures)
# THREE CASES:
# 1. TWO UTRS, ONE AT EITHER END - MOST COMMON
# 2. NO UTRS
# 3. A UTR AT ONE END BUT NOT THE OTHER - LABEL NEEDS TO KNOW WHICH END
# 
# NOW ALSO WORKS FOR MULTI-EXON GENES AND ALL EXAMPLES ON NEGATIVE STRAND.
  if ( length(disjoin(x)) == 1 ) {
    x = x
  } else if ( length(disjoin(x)) == 3 ) {
    x = unique(c(x, disjoin(x)))
    if ( unique(strand(x)) == "+" ) {
      names(x) = c("exon", "CDS", "5pUTR", "3pUTR")
    } else {
      names(x) = c("exon", "CDS", "3pUTR", "5pUTR")
    }
  } else if ( length(disjoin(x)) == 2 ) {
    if ( is.na(precede(setdiff(disjoin(x), x["CDS"]), x["CDS"])) ) {
      x = unique(c(x, disjoin(x)))
      names(x) = c("exon", "CDS", "3pUTR")
    } else {
      x = unique(c(x, disjoin(x)))
      names(x) = c("exon", "CDS", "5pUTR")
    }
  } else if ( length(disjoin(x)) > 3 ) {
    x = c(x, setdiff(disjoin(x), x[names(x) == "CDS"]))
    if ( length(setdiff(disjoin(x), x[names(x) == "CDS"])) == 2) {
      if ( unique(strand(x)) == "+" ) {
        names(x) = c(names(x)[1:(length(names(x))-2)], "5pUTR", "3pUTR")
      } else {
        names(x) = c(names(x)[1:(length(names(x))-2)], "3pUTR", "5pUTR")
      }
    } else {
      if ( is.na(precede(setdiff(disjoin(x), x[names(x) == "CDS"]), x[names(x) == "CDS"])) ) {
        names(x) = c(names(x)[1:(length(names(x))-1)], "3pUTR")
      } else {
        names(x) = c(names(x)[1:(length(names(x))-1)], "5pUTR")
      }
    }
    
  }
  return(x)
}
