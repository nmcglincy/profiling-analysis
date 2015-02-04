annotateSigChanges = function(x) {
  require("org.Sc.sgd.db")
  gene.name.l = as.list(org.Sc.sgdGENENAME)
  x$gene.symbol = unlist(gene.name.l[x$sys.id])
  gene.descriptions = as.list(org.Sc.sgdDESCRIPTION)
  x$gene.description = unlist(gene.descriptions[x$sys.id])
  x
}
# 
annotateYeastDF = function(x, y) {
	# where x is the DF and y is the column with the ORF ids
  require("org.Sc.sgd.db")
  gene.name.l = as.list(org.Sc.sgdGENENAME)
  x$gene.symbol = unlist(gene.name.l[unlist(x[, y])])
  gene.descriptions = as.list(org.Sc.sgdDESCRIPTION)
  x$gene.description = unlist(gene.descriptions[unlist(x[, y])])
  x
}