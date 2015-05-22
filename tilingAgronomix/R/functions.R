
###
### getPointscoords
plotDataGcAvg <- function(genename, exGainDelt = 1, inLossDelt = 1, exMinLength = 2.5, CONTROL = "CONTROL", exper = "brm", allGC = NULL, k=3) {
  tres1    <- getDataForGene(genename=genename, normMicData=matr1$rpb3_loess, type, 
                             positions, genes, geneNames, gffs, names, locsPos=locsPos, locsStrand=locsStrand, allGC=allGC) 
  if (is.null(tres1)) {
    return(NULL)
  }
  if (length(tres1$probesInAll) < 18) {
    return(NULL)
  }
  ktoryStrand <- (names(head(sort(table(tres1$exones[,7]), decreasing=TRUE),1)) == "+") + 0
  n.strand = length(cels)
  indyStr <- tres1$tsgen$strand == ktoryStrand
  if (sum(indyStr) < 18) {
    ktoryStrand <- 1 - ktoryStrand
    indyStr <- tres1$tsgen$strand == ktoryStrand
    if (sum(indyStr) < 18) {
      return(NULL)
    }
  }
  
  probsInAllOrAny =  (tres1$probesInAll[indyStr] | tres1$probesInAny[indyStr])[1:(sum(indyStr)/n.strand)]
  indForm <- cumsum(abs(c(0,diff(probsInAllOrAny))))
  
  sum1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), sum, na.rm=T)))
  mean1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), mean, na.rm=T)))
  length1 <- (sum1[,3]/mean1[,3])/3
  if (nrow(mean1) < 8) {
    return(NULL)
  }
  
  #
  # filters on WT
  exon1 =  tapply(probsInAllOrAny, indForm, any)
  tofilter <- data.frame(mean1[mean1[,2] == CONTROL,], 
                         exon =  exon1, 
                         length1 = length1[1:length(exon1)])
  tofilter$lengthT <- tofilter$length1 >= exMinLength
  tofilter$diffNei <- tofilter[,3] - pmax(tofilter[c(2:(nrow(tofilter)),nrow(tofilter)-1),3] , tofilter[c(2,1:(nrow(tofilter)-1)),3])
  tofilter$diffNeiT <- tofilter$diffNei >= exGainDelt
  tofilter$diffNieT <- tofilter$diffNei <= inLossDelt
  tofilter$exper <- mean1[mean1[,2] == exper,3]
  tofilter$experdiff <- mean1[mean1[,2] == exper,3] - mean1[mean1[,2] == CONTROL,3]
  
  # maska dla exonow
  #  maska <- rep(tofilter$diffNeiT & tofilter$lengthT & tofilter$exon, round(length1[1:length(exon1)]))
  maska <- rep(tofilter$exon, round(length1[1:length(exon1)]))
  frag <- rep(1:length(exon1), round(length1[1:length(exon1)]))
  tmp <- tres1$tsgen[indyStr, ]
  tdat <- data.frame(y = c(tmp[tmp[,4] == CONTROL,1], tmp[tmp[,4] == exper,1]),
                     probe = c(tmp[tmp[,4] == CONTROL,2], tmp[tmp[,4] == exper,2]),
                     frag = rep(frag, times=6),
                     exper = rep(c("control","experiment"), each=sum(tmp[,4] == exper))
  )
  
  gctmp <- tres1$tsgen$gc[indyStr]
  gc <- gctmp[seq(1,length(gctmp), 12)]
  tm <- unclass(by(tdat$y, list(tdat$probe, tdat$exper), mean, na.rm=TRUE))
  tm <- data.frame(tm, exon = rep(tofilter$exon, round(length1[1:length(exon1)])),
                   length = rep(tofilter$length1, round(length1[1:length(exon1)])),
                   diffNei = rep(tofilter$diffNei, round(length1[1:length(exon1)])),
                   exper = rep(tofilter$exper, round(length1[1:length(exon1)])),
                   experdiff = rep(tofilter$experdiff, round(length1[1:length(exon1)])),
                   gc = gc)
  tm
}

###
plotDtotPlot <- function(genename, exGainDelt = 1, inLossDelt = 1, exMinLength = 2.5, CONTROL = "CONTROL", exper = "brm", allGC = NULL, k=3) {
  tres1    <- getDataForGene(genename=genename, normMicData=matr1$rpb3_loess, type, 
                             positions, genes, geneNames, gffs, names, locsPos=locsPos, locsStrand=locsStrand, allGC=allGC) 
  if (is.null(tres1)) {
    return(NULL)
  }
  if (length(tres1$probesInAll) < 18) {
    return(NULL)
  }
  ktoryStrand <- (names(head(sort(table(tres1$exones[,7]), decreasing=TRUE),1)) == "+") + 0
  n.strand = length(cels)
  indyStr <- tres1$tsgen$strand == ktoryStrand
  if (sum(indyStr) < 18) {
    ktoryStrand <- 1 - ktoryStrand
    indyStr <- tres1$tsgen$strand == ktoryStrand
    if (sum(indyStr) < 18) {
      return(NULL)
    }
  }
  
  probsInAllOrAny =  (tres1$probesInAll[indyStr] | tres1$probesInAny[indyStr])[1:(sum(indyStr)/n.strand)]
  indForm <- cumsum(abs(c(0,diff(probsInAllOrAny))))
  
  sum1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), sum, na.rm=T)))
  mean1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), mean, na.rm=T)))
  length1 <- (sum1[,3]/mean1[,3])/3
  if (nrow(mean1) < 8) {
    return(NULL)
  }
  
  #
  # filters on WT
  exon1 =  tapply(probsInAllOrAny, indForm, any)
  tofilter <- data.frame(mean1[mean1[,2] == CONTROL,], 
                         exon =  exon1, 
                         length1 = length1[1:length(exon1)])
  tofilter$lengthT <- tofilter$length1 >= exMinLength
  tofilter$diffNei <- tofilter[,3] - pmax(tofilter[c(2:(nrow(tofilter)),nrow(tofilter)-1),3] , tofilter[c(2,1:(nrow(tofilter)-1)),3])
  tofilter$diffNeiT <- tofilter$diffNei >= exGainDelt
  tofilter$diffNieT <- tofilter$diffNei <= inLossDelt
  tofilter$exper <- mean1[mean1[,2] == exper,3]
  tofilter$experdiff <- mean1[mean1[,2] == exper,3] - mean1[mean1[,2] == CONTROL,3]
  
  # maska dla exonow
  #  maska <- rep(tofilter$diffNeiT & tofilter$lengthT & tofilter$exon, round(length1[1:length(exon1)]))
  maska <- rep(tofilter$exon, round(length1[1:length(exon1)]))
  frag <- rep(1:length(exon1), round(length1[1:length(exon1)]))
  tmp <- tres1$tsgen[indyStr, ]
  tdat <- data.frame(y = c(tmp[tmp[,4] == CONTROL,1][maska], tmp[tmp[,4] == exper,1][maska]),
                     probe = c(tmp[tmp[,4] == CONTROL,2][maska], tmp[tmp[,4] == exper,2][maska]),
                     frag = rep(frag[maska], times=6),
                     exper = rep(c("control","experiment"), each=sum(maska)*3))
  
  tm<-unclass(by(tdat$y, list(tdat$probe, tdat$exper), mean))
  #   plot(tm[,1], tm[,2],xlim=c(5,13), ylim=c(5,13), xlab=colnames(tm)[1], ylab=colnames(tm)[2], main=genename)
  #   abline(0,1)
  # 
  #   plot(tm[,1], tm[,2]-tm[,1],xlim=c(5,13), ylim=c(-3,3), xlab=colnames(tm)[1], ylab="diff", main=genename)
  #   abline(h=0)
  #   
  library(zoo)
  tm2 <- apply(tm, 2, rollmean,k=k)
  plot(tm2[,1], tm2[,2]-tm2[,1],xlim=c(5,13), ylim=c(-3,3), xlab=colnames(tm)[1], ylab="diff", main=genename, type="b",pch=19,cex=0.7)
  abline(h=0)  
}

###
### p -values
calculatePValuesExonIntron <- function(genename, exGainDelt = 1, inLossDelt = 1, exMinLength = 2.5, CONTROL = "CONTROL", exper = "brm", allGC = NULL) {
  tres1    <- getDataForGene(genename=genename, normMicData=matr1$rpb3_loess, type, 
                             positions, genes, geneNames, gffs, names, locsPos=locsPos, locsStrand=locsStrand, allGC=allGC) 
  if (is.null(tres1)) {
    return(NULL)
  }
  if (length(tres1$probesInAll) < 18) {
    return(NULL)
  }
  ktoryStrand <- (names(head(sort(table(tres1$exones[,7]), decreasing=TRUE),1)) == "+") + 0
  n.strand = length(cels)
  indyStr <- tres1$tsgen$strand == ktoryStrand
  if (sum(indyStr) < 18) {
    ktoryStrand <- 1 - ktoryStrand
    indyStr <- tres1$tsgen$strand == ktoryStrand
    if (sum(indyStr) < 18) {
      return(NULL)
    }
  }
  
  probsInAllOrAny =  (tres1$probesInAll[indyStr] | tres1$probesInAny[indyStr])[1:(sum(indyStr)/n.strand)]
  indForm <- cumsum(abs(c(0,diff(probsInAllOrAny))))
  
  sum1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), sum, na.rm=T)))
  mean1 <- as.data.frame(as.table(by(tres1$tsgen[indyStr, 1], list( part = rep(indForm, n.strand), mutant = tres1$tsgen[indyStr, "mutant"]), mean, na.rm=T)))
  length1 <- (sum1[,3]/mean1[,3])/3
  if (nrow(mean1) < 8) {
    return(NULL)
  }
  
  #
  # filters on WT
  exon1 =  tapply(probsInAllOrAny, indForm, any)
  tofilter <- data.frame(mean1[mean1[,2] == CONTROL,], 
                         exon =  exon1, 
                         length1 = length1[1:length(exon1)])
  tofilter$lengthT <- tofilter$length1 >= exMinLength
  tofilter$diffNei <- tofilter[,3] - pmax(tofilter[c(2:(nrow(tofilter)),nrow(tofilter)-1),3] , tofilter[c(2,1:(nrow(tofilter)-1)),3])
  tofilter$diffNeiT <- tofilter$diffNei >= exGainDelt
  tofilter$diffNieT <- tofilter$diffNei <= inLossDelt
  tofilter$exper <- mean1[mean1[,2] == exper,3]
  tofilter$experdiff <- mean1[mean1[,2] == exper,3] - mean1[mean1[,2] == CONTROL,3]
  
  # maska dla exonow
  maska <- rep(tofilter$diffNeiT & tofilter$lengthT & tofilter$exon, round(length1[1:length(exon1)]))
  frag <- rep(1:length(exon1), round(length1[1:length(exon1)]))
  tmp <- tres1$tsgen[indyStr, ]
  tdat <- data.frame(y = c(tmp[tmp[,4] == CONTROL,1][maska], tmp[tmp[,4] == exper,1][maska]),
                     probe = c(tmp[tmp[,4] == CONTROL,2][maska], tmp[tmp[,4] == exper,2][maska]),
                     frag = rep(frag[maska], times=6),
                     exper = rep(c("control","experiment"), each=sum(maska)*3))
  
  tm<-unclass(by(tdat$y, list(tdat$probe, tdat$exper), mean))
  plot(tm[,1], tm[,2])
  
  
  exon.p.vals <- try(anova(model <- lm(y ~ factor(frag) * exper + factor(probe), data=tdat))[c(2,4),5], silent=T)
  if (class(exon.p.vals) == "try-error") {
    exon.p.vals = NULL
    mainexon = 0
  } else {
    mainexon= model$coefficients["experexperiment"]
  }
  
  # maska dla intronow
  maska <- rep(tofilter$diffNieT & tofilter$lengthT & !tofilter$exon, round(length1[1:length(exon1)]))
  tdat2 <- data.frame(y = c(tmp[tmp[,4] == CONTROL,1][maska], tmp[tmp[,4] == exper,1][maska]),
                      probe = c(tmp[tmp[,4] == CONTROL,2][maska], tmp[tmp[,4] == exper,2][maska]),
                      frag = rep(frag[maska], times=6),
                      exper = rep(c("control","experiment"), each=sum(maska)*3))
  intron.p.vals <- try(anova(model2 <- lm(y ~ factor(frag) * exper + factor(probe), data=tdat2))[c(2,4),5], silent=T)
  if (class(intron.p.vals) == "try-error") {
    intron.p.vals = NULL
    mainintron = 0
  } else {
    mainintron= model2$coefficients["experexperiment"]
  }
  list(exon.p.vals = exon.p.vals, intron.p.vals = intron.p.vals,
       diff.exo = diff(range(tofilter[tofilter$diffNeiT & tofilter$lengthT & tofilter$exon,"experdiff"])),
       diff.int = diff(range(tofilter[tofilter$diffNieT & tofilter$lengthT & !tofilter$exon,"experdiff"])),
       mainexon=mainexon, mainintron=mainintron)
}


getAndPlotGene <- function(genename, which = -3,genename2=genename, ...) {
  tres1 <- getDataForGene(genename=genename, normMicData=matr1$rpb3_loess, type, 
                          positions, genes, geneNames, gffs, names, locsPos=locsPos, locsStrand=locsStrand) 
  if(length(tres1) ==8) {
    ktoryStrand <- (names(head(sort(table(tres1$exones[,7]), decreasing=TRUE),1)) == "+") + 0
    indyStr <- tres1$tsgen$strand == ktoryStrand
    if (sum(indyStr[1:length(tres1$gidds)]) < 12) indyStr = !indyStr
    if (sum(indyStr[1:length(tres1$gidds)]) < 12) return(0)
    
    #     if (which == -3)
    #       indyStr[(2*length(indyStr)/3 + 1):(length(indyStr))] = FALSE
    #     if (which == -2)
    #       indyStr[(length(indyStr)/3 + 1):(2*length(indyStr)/3)] = FALSE
    layout(matrix(1:2,2,1),heights=c(2,1))
    par(mar=c(0,0,0,0), oma=c(2,3,4,2))
    plotGene(tres1$tsgen[indyStr,], 
             tres1$exones, 
             tres1$gidds[indyStr[1:length(tres1$gidds)]], 
             tres1$probesInAll[indyStr], 
             tres1$probesInAny[indyStr], 
             type, 
             tres1$isoformsCode[indyStr], genename=genename2,whichaxis=3, ...) 
    
    plotGeneDiff(tres1$tsgen[indyStr,], tres1$exones, tres1$gidds[indyStr[1:length(tres1$gidds)]], 
                 tres1$probesInAll[indyStr], tres1$probesInAny[indyStr], 
                 type, tres1$isoformsCode[indyStr], genename="", whichaxis=1, plotisoforms=FALSE, ...) 
  }
}


plotGene <- function(tsgen, exones, gidds, probesInAll, probesInAny, type, isoformsCode, minNoProbes = 12, genename="",whichaxis=1, plotAll=F, ncsel = NULL) {
  avgs <- unclass(by(tsgen$y, list(tsgen$probe, tsgen$mutant), mean, na.rm=T))
  filtr <- apply(is.na(avgs), 2, mean) < 0.5
  avgs <- avgs[,filtr, drop=F]
  cols <- brewer.pal(2*ncol(avgs), "Paired")
  
  matplot(c(positions[gidds]-10,min(positions[gidds])+2000), rbind(avgs,10), type="n",las=1,lty=1,xlab="",ylab="",
          ylim=c(4,14),xaxt="n")
  mtext(genename, 3, line=3)
  axis(whichaxis)
  #
  # the background
  rect(positions[gidds][probesInAll], 3, positions[gidds][probesInAll]+25, 20, 
       col="#AAAAAA11", border="#AAAAAA11")
  if (sum(xor(probesInAny,probesInAll)) > 0)
    rect(positions[gidds][xor(probesInAny,probesInAll)], 3, positions[gidds][xor(probesInAny,probesInAll)]+25, 20, 
         col="#77777733", border="#77777733")
  #
  # probes
  if (is.null(ncsel)) 
    ncsel <- 1:ifelse(plotAll, ncol(avgs), min(2,ncol(avgs)))
  ncsel <- intersect(ncsel, which(filtr))
  if (length(ncsel) < 1) return(0)
  for (nc in ncsel) { #ncol(avgs)
    lines(positions[gidds]+12, lowess(positions[gidds]+12, avgs[,nc], min(0.5,7/nrow(avgs)))$y, col=cols[nc*2-1])
    
    lines(rep(positions[gidds], each=2) + c(0,25), rep(avgs[,nc], each=2), col=cols[nc*2-1], lwd=0.5, lty="29")
    rect(positions[gidds], avgs[,nc], positions[gidds]+25, avgs[,nc], col=cols[nc*2], border=cols[nc*2], lwd=3)
  }
  legend("topleft", unique(type)[ncsel], lwd=3, col=cols[(1:(length(cols)/2))*2][ncsel], bty="n", ncol=length(cols)/2, cex=0.8)
  llab <- c("exon","five_prime_UTR","three_prime_UTR")
  lcols <- c("#AAAAAA77","#0000AA77","#AA000077")
  rect(exones[,4], 
       3.7+as.numeric(factor(exones[,9]))/5,
       exones[,5], 
       3.8+as.numeric(factor(exones[,9]))/5,
       col=lcols[as.numeric(factor(exones[,3], levels=llab))], 
       border=lcols[as.numeric(factor(exones[,3], levels=llab))])
  
}

plotGeneDiff <- function(tsgen, exones, gidds, probesInAll, probesInAny, type, isoformsCode, minNoProbes = 12, genename="", whichaxis=1, plotisoforms=TRUE, plotAll=FALSE, ncsel = NULL) {
  avgs <- unclass(by(tsgen$y, list(tsgen$probe, tsgen$mutant), mean))
  filtr <- apply(is.na(avgs), 2, mean) < 0.5
  avgs <- avgs[,filtr, drop=F]
  cols <- brewer.pal(2*ncol(avgs), "Paired")
  
  matplot(c(positions[gidds]-10,min(positions[gidds])+2000), rbind(avgs,10), type="n",las=1,lty=1,xlab="",ylab="",
          main=genename, ylim=c(-3.5,3.5),xaxt="n",yaxt="n")
  axis(whichaxis)
  axis(4,las=1)
  abline(h=0)
  abline(h=-4:4, lty=3, col="grey")
  #
  # the background
  rect(positions[gidds][probesInAll], -20, positions[gidds][probesInAll]+25, 20, 
       col="#AAAAAA11", border="#AAAAAA11")
  if (sum(xor(probesInAny,probesInAll)) > 0)
    rect(positions[gidds][xor(probesInAny,probesInAll)], -20, positions[gidds][xor(probesInAny,probesInAll)]+25, 20, 
         col="#77777733", border="#77777733")
  #
  # probes
  if (is.null(ncsel)) 
    ncsel <- 1:ifelse(plotAll, ncol(avgs), min(2,ncol(avgs)))
  ncsel <- intersect(ncsel, which(filtr))
  if (length(ncsel) < 2) return(0)
  for (nc in ncsel[-1]) { #ncol(avgs)
    lines(positions[gidds]+12, lowess(positions[gidds]+12, avgs[,nc]-avgs[,ncsel[1]], min(0.5,7/nrow(avgs)))$y, col=cols[nc*2-1])
    
    lines(rep(positions[gidds], each=2) + c(0,25), rep(avgs[,nc]-avgs[,ncsel[1]], each=2), col=cols[nc*2-1], lwd=0.5, lty="29")
    rect(positions[gidds], avgs[,nc]-avgs[,ncsel[1]], positions[gidds]+25, avgs[,nc]-avgs[,ncsel[1]], col=cols[nc*2], border=cols[nc*2], lwd=3)
  }
  legend("topright", paste(unique(type)[ncsel][-1],"-",unique(type)[ncsel[1]]), lwd=3, col=cols[(2:(length(cols)/2))*2], bty="n", ncol=length(cols)/2-1, cex=0.8)
  llab <- c("exon","five_prime_UTR","three_prime_UTR")
  lcols <- c("#AAAAAA77","#0000AA77","#AA000077")
  if(plotisoforms) {
    rect(exones[,4], 
         -2.1+as.numeric(factor(exones[,9]))/5,
         exones[,5], 
         -2+as.numeric(factor(exones[,9]))/5,
         col=lcols[as.numeric(factor(exones[,3], levels=llab))], 
         border=lcols[as.numeric(factor(exones[,3], levels=llab))])
  }
  
}






loadTilingArrays <- function (cels, names, type, bpmap, bpremap, method="loess", AGRONOMIX=FALSE) {
  bpmapChr1 <- readBpmap(bpmap)
  chrs <- c("chrC", "chr1", "chr2", "chr3", "chr4", "chr5", "chrM")
  if(AGRONOMIX) {
    # AGRONOMIX
    bpmapChr2 <- bpmapChr1[c(68,63:67,69)]
  }else {
    bpmapChr2 <- bpmapChr1[61:67]
  }
  
  locs <- read.table(bpremap, sep="\t")
  for (kt in 2:6) {
    tmp <- locs[locs[,1]==chrs[kt],]
    bpmapChr2[[kt]]$seqInfo$numberOfHits = nrow(tmp)
    bpmapChr2[[kt]]$pmx = as.integer(tmp[,5])
    bpmapChr2[[kt]]$pmy = as.integer(tmp[,6])
    bpmapChr2[[kt]]$probeseq = as.character(tmp[,4])
    bpmapChr2[[kt]]$startpos = as.integer(tmp[,2])
    bpmapChr2[[kt]]$strand =  as.integer(tmp[,3])
  }
  
  rpb3Chr1 <- readCelFile(bpmapChr2, cels, names, type, featureData=T, log.it=T)
  #  rpb3Chr1F <- readCelFile(bpmapChr2, cels, names, type, featureData=T, log.it=F)
  rpb3_loess <- normalize.Probes(rpb3Chr1, method=method)
  #  rpb3_loess = rpb3Chr1
  
  positions <- featureData(rpb3_loess)$pos
  if(AGRONOMIX) {
    # AGRONOMIX
    chromosomes <-  as.numeric(substr(featureData(rpb3_loess)$chr,13,13)) 
  }else {
    chromosomes <-  as.numeric(substr(featureData(rpb3_loess)$chr,14,14))
  }
  chromosomes[is.na(chromosomes)] <- 0
  
  list(rpb3_loess=rpb3_loess, 
       positions=positions, chromosomes=chromosomes, locs = locs[,1:3])
}


getDataForGene <- function(genename, normMicData, type, positions, genes, geneNames, gffs, names=names, probelength=25, probemargin=5, locsPos, locsStrand, allGC = NULL) {
  geneid <- head(which(geneNames == genename),1)
  if (length(geneid) == 0) return(NULL)
  #
  # probes in the gene region
  start <- genes[geneid, 4]
  stop  <- genes[geneid, 5]  
  chr   <- as.character(genes[geneid,1])
  gidds <- which((positions >= start - probemargin) & (positions <= stop - probelength+probemargin) & (chromosomes == as.numeric(substr(chr,4,4))))
  
  if (length(gidds) > 0) {
    gidds2 <- which((locsPos[[tolower(chr)]] >= start - probemargin) & (locsPos[[tolower(chr)]] <= stop - probelength+probemargin))
    if (length(gidds2) == length(gidds)) {
      strand <- locsStrand[[tolower(chr)]][gidds2]      
    } else {
      strand <- rep(2, length(gidds))
    }
    tx <- exprs(normMicData)[gidds,,drop=F]
    y  <- as.vector(tx)
    probe  <- rep(1:nrow(tx), ncol(tx))
    probePos <- rep(positions[gidds], ncol(tx))
    probeStrand <- rep(strand, ncol(tx))
    mutant <- factor(rep(type, each=nrow(tx)), levels=unique(type))
    plate  <- rep(names, each=nrow(tx))
    
    gcAlong <- 0
    if (!is.null(allGC)) {
      gcAlong <- rep(allGC[gidds], each=ncol(tx))
    }
    
    tsgen  <- data.frame(y, probe, probePos, mutant, plate, 
                         gm =  strsplit(as.character(genes[geneid,9]), split="Name=")[[1]][2],
                         strand = probeStrand, gc = gcAlong)
  } else {
    return(NULL)
  }
  
  #
  # exones in the gene region
  tmp <- gffs[gffs[,4] >= start & gffs[,5] <= stop,,drop=F]
  tmp <- tmp[grep(tmp[,9], pattern=geneNames[geneid]),,drop=F]
  exones <- tmp[tmp[,3] %in% c("five_prime_UTR", "exon", "three_prime_UTR"),,drop=F]
  if (nrow(exones) == 0 )
    return(NULL)
  # gene models
  lev <- levels(factor(exones[,9]))
  # coding of gene models 
  coding <- sapply(lev, function(le) {
    etmp <- exones[exones[,9] == le,]
    apply(apply(etmp[,4:5],1,function(coord) 
      tsgen[,3] >= (coord[1] - probemargin) & (tsgen[,3] <= coord[2] - probelength + probemargin)
    ), 1, any)+0
  })
  colnames(coding) <- paste("F",substr(colnames(coding),18,20),sep="")
  
  #
  # probes common in all isoforms
  probesInAll  <- apply(coding==1, 1, any)
  probesInAny  <- apply(coding==1, 1, all)
  isoformsCode <- apply(coding, 1, function(x) sum(x*2^(1:length(x))))
  
  return(list(genename = genename, geneid = geneid,gidds=gidds,
              tsgen = tsgen, 
              exones = exones, 
              isoformsCode = isoformsCode, probesInAll=probesInAll, probesInAny=probesInAny,
              tx = tx))
}


# calculate change only in utrs
calculatePvaluesForUTR <- function (tsgen, exones, probelength=25, probemargin=5) {
  utrs <- c("five_prime_UTR", "three_prime_UTR")
  utrPvals <- lapply(utrs, function(reg) {
    primePvals <- NULL
    prime  <- exones[exones[,3] == reg & as.numeric(factor(exones[,9])) == 1, , drop=F]
    if(nrow(prime) >0) {
      wind <- apply(sapply(1:nrow(prime), function(nr) {
        (tsgen$probePos > prime[nr,4] - probemargin) & (tsgen$probePos < prime[nr,5] - probelength+probemargin) 
      }), 1, any)
      tmp <- tsgen[wind ,]
      if (nrow(tmp) > 0) {
        primePvals  <- 
          if (length(unique(tmp$probe)) == 1) 
            summary(lm(y ~ mutant, data=tmp))$coef[seq_along(unique(type)),c(1,4)]
        else
          summary(lm(y ~ mutant+factor(probe), data=tmp))$coef[seq_along(unique(type)),c(1,4)]
      }
    }
    primePvals
  })
  names(utrPvals) <- utrs
  utrPvals
}

calculatePvaluesForDiffs <- function (tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) {
  n.samples = nlevels(tsgen$plate)
  #
  # p-values that differences between isoforms are significant
  if (sum(probesInAny != probesInAll) == n.samples) {   # one probe only in diff
    form <- (y - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
               factor(mutant))
  } else {                               # many probes in diff
    labels <- rep(c(0,cumsum(diff((probesInAll+probesInAny)[1:(length(isoformsCode)/n.samples)]) !=0)), n.samples)
    number.of.exones.in.diff <- length(unique(labels[probesInAll != probesInAny]))
    if (number.of.exones.in.diff == 1) { # only one exon in diff
      form  <- (y - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
                  factor(probe)+factor(mutant)-1)
    } else {                            # many exons in diff
      form  <- (y - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
                  factor(probe)+factor(labels):factor(mutant)-1)
    }
  }
  
  pcoef  <- summary(lm(form, data=tsgen, subset= probesInAny != probesInAll))$coef
  pexones <- tail(anova(lm(form, data=tsgen, subset= probesInAny != probesInAll))[,5],2)[1] # one before last p-value
  list(pcoef=pcoef, pexones=pexones)
}

calculatePvaluesForIntrons2 <- function (tsgen, meansEffect, probesInAll, probesInAny, isoformsCode, mutantID="splicesom") {
  n.samples = nlevels(tsgen$plate)
  tmp <- !(probesInAll | probesInAny)[1:(length(probesInAll)/n.samples)]
  probesLabels <- rep(cumsum(!tmp), n.samples)
  probesIntrone <- rep(tmp, n.samples)
  
  number.of.exones.in.diff <- length(unique(probesLabels[probesIntrone]))
  if (number.of.exones.in.diff < 1) return(-2)
  
  probesLabels2 <- probesLabels
  probesLabels2[!tmp] <- -1
  #  probesLabels2[!(probesLabels %in% unique(probesLabels[probesIntrone]))] = -1
  form  <- (y ~ factor(probesLabels2)*factor(mutant == mutantID))
  summary(lm(form, data=tsgen))$coef[,c(1,2,4)]
}

calculatePvaluesForIntrons <- function (tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) {
  n.samples = nlevels(tsgen$plate)
  
  # is it intron
  tmp <- !(probesInAll | probesInAny)[1:(length(probesInAll)/n.samples)]
  probesLabels <- rep(cumsum(!tmp), n.samples)
  probesIntrone <- rep(tmp, n.samples)
  #
  # p-values that differences between isoforms are significant
  if (sum(probesIntrone) == n.samples) {   # one probe only in diff
    form <- (y ~ factor(mutant))
  } else {                               # many probes in diff
    number.of.exones.in.diff <- length(unique(probesLabels[probesIntrone]))
    if (number.of.exones.in.diff == 1) { # only one exon in diff
      form  <- (y ~ factor(probe)+factor(mutant)-1)
    } else {                            # many exons in diff
      form  <- (y ~ factor(probe)+factor(probesLabels):factor(mutant)-1)
    }
  }
  
  pcoef  <- summary(lm(form, data=tsgen, subset= probesIntrone))$coef
  pexones <- tail(anova(lm(form, data=tsgen, subset= probesIntrone))[,5],2)[1] # one before last p-value
  
  list(pcoef=pcoef, pexones=pexones)
}
calculatePvaluesForExones <- function (tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) {
  n.samples = nlevels(tsgen$plate)
  
  # is it exon
  tmp <- (probesInAll | probesInAny)[1:(length(probesInAll)/n.samples)]
  probesLabels <- rep(cumsum(!tmp), n.samples)
  probesIntrone <- rep(tmp, n.samples)
  #
  # p-values that differences between isoforms are significant
  if (sum(probesIntrone) == n.samples) {   # one probe only in diff
    form <- (y  - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
               factor(mutant))
  } else {                               # many probes in diff
    number.of.exones.in.diff <- length(unique(probesLabels[probesIntrone]))
    if (number.of.exones.in.diff == 1) { # only one exon in diff
      form  <- (y  - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
                  factor(probe)+factor(mutant)-1)
    } else {                            # many exons in diff
      form  <- (y  - meansEffect[paste("mutant",mutant,sep=""), 1] ~ 
                  factor(probe)+factor(probesLabels):factor(mutant)-1)
    }
  }
  
  pcoef  <- summary(lm(form, data=tsgen, subset= probesIntrone))$coef
  pexones <- tail(anova(lm(form, data=tsgen, subset= probesIntrone))[,5],2)[1] # one before last p-value
  
  list(pcoef=pcoef, pexones=pexones)
}

getMaxEffect <- function(tmp) {
  tmp2 <- tmp[-grep(rownames(tmp),pattern="\\(probe\\)"),,drop=F]
  tmp2 <- tmp2[-grep(rownames(tmp2),pattern="Intercept"),,drop=F]
  tmp2[which.max(abs(tmp2[,1])),1]
}

# assume that there are only two groups to be compared
# meansEffect - main differences between genes calculated on common parts
# 
#
calculatePvaluesForGene2groups <- function(tsgen, exones, probesInAll, probesInAny, type, 
                                           isoformsCode, minNoProbes = 12,
                                           whichPV = c("main","intrones","intrones2","exones","diffsInExones","utr"),
                                           mutantID="splicesom") {
  if (sum(probesInAll) < minNoProbes) return(NULL)
  
  n.samples = nlevels(tsgen$plate)
  number.of.probes.in.diff <- sum(probesInAny != probesInAll)/n.samples
  number.of.forms <- length(levels(factor(exones[,9])))
  
  results <- NULL
  if ("main" %in% whichPV) {
    meansPvalues <- summary(lm(y ~ mutant+factor(probe), data=tsgen, subset= probesInAll))$coef[seq_along(unique(type)),c(1,4)]
    meansEffect  <- summary(lm(y ~ mutant+factor(probe)-1, data=tsgen, subset= probesInAll))$coef[seq_along(unique(type)),c(1,4)]
    
    results$main <- list(meansPvalues=meansPvalues, meansEffect=meansEffect)
  }
  
  if ("intrones" %in% whichPV) {
    pvaluesForIntron <- NULL
    maxEffect <- NULL
    if (sum(!(probesInAll | probesInAny)) > 0) {
      pvaluesForIntron <- calculatePvaluesForIntrons(tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) 
      maxEffect <- getMaxEffect(pvaluesForIntron[[1]])
    } 
    results$intrones <- list(pvaluesForIntron=pvaluesForIntron, intronMaxEffect=maxEffect)
  }
  
  if ("intrones2" %in% whichPV) {
    pvaluesForIntron2 <- NULL
    if (sum(!(probesInAll | probesInAny)) > 0) {
      pvaluesForIntron2 <- calculatePvaluesForIntrons2(tsgen, meansEffect, probesInAll, probesInAny, isoformsCode, mutantID=mutantID) 
    } 
    results$intrones2 <- list(pvaluesForIntron2=pvaluesForIntron2)
  }
  
  if ("exones" %in% whichPV) {
    pvaluesForExon <- NULL
    maxEffect <- NULL
    if (sum(!(probesInAll | probesInAny)) > 0) {
      pvaluesForExon <- calculatePvaluesForExones(tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) 
      maxEffect <- getMaxEffect(pvaluesForExon[[1]])
    } 
    results$exones <- list(pvaluesForExon=pvaluesForExon, exonMaxEffect=maxEffect)
  }
  
  if ("diffsInExones" %in% whichPV) {
    pvaluesForDiffs <- NULL
    maxEffect <- NULL
    if (number.of.forms>1 & number.of.probes.in.diff > 0) {
      pvaluesForDiffs <- calculatePvaluesForDiffs(tsgen, meansEffect, probesInAll, probesInAny, isoformsCode) 
      maxEffect <- getMaxEffect(pvaluesForDiffs[[1]])
    } 
    results$diffsInExones <- list(pvaluesForDiffs=pvaluesForDiffs, diffsMaxEffect=maxEffect)
  }
  
  if ("utr" %in% whichPV) {
    pvaluesForUTR <- NULL
    if (sum(!(probesInAll | probesInAny)) > 0) {
      pvaluesForUTR <- calculatePvaluesForUTR(tsgen, exones)
    } 
    results$utr <- list(pvaluesForUTR=pvaluesForUTR)
  }
  
  results$length = sum(probesInAll)/n.samples
  results
}


getVals <- function(plist, fun, naval) {
  sapply(plist, function(x) {
    tmp <- eval(parse(text=fun))
    ifelse(is.null(tmp), naval, tmp)
  })
}

