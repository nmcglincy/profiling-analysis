UtrCreator = function(x) {
  require(GenomicRanges)
  require(GenomicFeatures)
# THREE CASES:
# 1. TWO UTRS, ONE AT EITHER END - MOST COMMON
# 2. NO UTRS
# 3. A UTR AT ONE END BUT NOT THE OTHER - LABEL NEEDS TO KNOW WHICH END
  if ( length(disjoin(x)) == 1 ) {
    x = x
  } else if ( length(disjoin(x)) == 3 ) {
    x = unique(c(x, disjoin(x)))
    names(x) = c("exon", "CDS", "5pUTR", "3pUTR")
  } else if ( length(disjoin(x)) == 2 ) {
    if ( follow(disjoin(x)[2], disjoin(x)[1]) == 1 ) {
      x = unique(c(x, disjoin(x)))
      names(x) = c("exon", "CDS", "3pUTR")
    } else {
      x = unique(c(x, disjoin(x)))
      names(x) = c("exon", "CDS", "5pUTR")
    }
  }
  return(x)
}
