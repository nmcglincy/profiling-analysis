UtrCreator = function(x) {
  require(GenomicRanges)
  require(GenomicFeatures)
  if ( length(disjoin(x)) == length(x)/2 ) {
    foo = x 
    mcols(foo) = NULL
    x = foo
  } else if ( length(disjoin(x)) == 3 ) {
    x = unique(c(x, disjoin(x), ignore.mcols = TRUE))
    if ( unique(strand(x)) == "+" ) {
      names(x) = c("exon", "CDS", "5pUTR", "3pUTR")
    } else {
      names(x) = c("exon", "CDS", "3pUTR", "5pUTR")
    }
  } else if ( length(disjoin(x)) == 2 ) {
    if ( is.na(precede(setdiff(disjoin(x), x["CDS"]), x["CDS"])) ) {
      x = unique(c(x, disjoin(x), ignore.mcols = TRUE))
      names(x) = c("exon", "CDS", "3pUTR")
    } else {
      x = unique(c(x, disjoin(x), ignore.mcols = TRUE))
      names(x) = c("exon", "CDS", "5pUTR")
    }
  } else if ( length(disjoin(x)) > 3 ) {
    x = c(x, setdiff(disjoin(x), x[names(x) == "CDS"]), ignore.mcols = TRUE)
    if (length(setdiff(disjoin(x), x[names(x) == "CDS"])) == 0) {
      x = x
    } else if ( length(setdiff(disjoin(x), x[names(x) == "CDS"])) == 2) {
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
