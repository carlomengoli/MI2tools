getAllProbesForGene <- function(genename, positions, chromosomes, genes, geneNames, probelength=25, probemargin=5) {
  geneid <- head(which(geneNames == genename),1)
  if (length(geneid) == 0) return(NULL)
  #
  # probes in the gene region
  start <- genes[geneid, 4]
  stop  <- genes[geneid, 5]  
  chr   <- as.character(genes[geneid,1])
  gidds <- which((positions >= start - probemargin) & (positions <= stop - probelength+probemargin) & (chromosomes == as.numeric(substr(chr,4,4))))
  gidds
}
