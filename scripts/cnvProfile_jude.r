library("DNAcopy", lib.loc="/home/cnv/CNV_dataprocessing/DNAcopy_1.26.0")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

remove.segment <- function( rsShort, rsSegnum, ratioData, sd.undo ) {

	appendLeft <- TRUE
	checkSdundo <- FALSE

	if (rsSegnum == 1) {
		appendLeft <- FALSE
	} else {
	if (rsSegnum == nrow(rsShort)) {
		appendLeft <- TRUE
	} else {
		rightIndex <- rsSegnum + 1
		leftIndex <- rsSegnum - 1
		
		if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
			appendLeft <- TRUE
		} else {
		if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
			appendLeft <- FALSE
		} else {
		if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) < abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
			appendLeft <- TRUE
			checkSdundo <- TRUE
		} else {
			appendLeft <- FALSE
			checkSdundo <- TRUE
		}}}
	}}
	
	appendIndex <- 99999999
	if (appendLeft) {
		appendIndex <- rsSegnum - 1
	} else {
		appendIndex <- rsSegnum + 1
	}
	
	tempShort <- rsShort
	newLocStart <- -1
	newLocEnd <- -1
	if (appendLeft) {
		tempShort[appendIndex, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
		tempShort[appendIndex, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
	} else {
		tempShort[appendIndex, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
		tempShort[appendIndex, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
	}
	
	tempShort[appendIndex, "num.mark"] <- tempShort[appendIndex, "num.mark"] + tempShort[rsSegnum, "num.mark"]
	tempShort[appendIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[appendIndex, "seg.start"]:tempShort[appendIndex, "seg.end"]], base=2))
	
	cat("append", tempShort[appendIndex, "chrom"], tempShort[appendIndex, "loc.start"], tempShort[appendIndex, "loc.end"], tempShort[appendIndex, "num.mark"], tempShort[appendIndex, "seg.mean"], tempShort[appendIndex, "seg.start"], tempShort[appendIndex, "seg.end"], "\n")
	
	tempShort <- tempShort[-rsSegnum, ]
	tempShort$segnum <- seq(1:nrow(tempShort))
	
	if (checkSdundo) {
		thisSd <- -1
		if (appendLeft) {
			leftIndex <- appendIndex
			rightIndex <- appendIndex + 1
		} else {
			leftIndex <- appendIndex - 2
			rightIndex <- appendIndex - 1
		}
		#thisSd <- sd(ratioData[tempShort$seg.start[leftIndex]:tempShort$seg.start[rightIndex], "lowratio"])
		thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
		
		if (abs(tempShort$seg.mean[leftIndex] - tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {

			cat("left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
			cat("right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

			##  remove breakpoint
			tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
			tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
			tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
			tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
			tempShort <- tempShort[-rightIndex, ]
			tempShort$segnum <- seq(1:nrow(tempShort))
		}
	}
	
	return(tempShort)
}


sdundo.all <- function (sdShort, ratioData, sd.undo) {

	tempShort <- sdShort
	thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)
	
	while ( TRUE ) {
		
		chrom <- tempShort$chrom
		chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])
	
		breakpoints <- which(chrom == chrom.shift)
		cat("sdundo.all intrachrom breakpoints", length(breakpoints), "\n")
		
		if (length(breakpoints) < 1) {
			break
		}
		
		breakpoints.shift <- breakpoints + 1
				
		undo.breakpoints <- breakpoints[which(abs(tempShort$seg.mean[breakpoints] - tempShort$seg.mean[breakpoints.shift]) < thisSd * sd.undo)]

		cat("sdundo.all undo breakpoints", length(undo.breakpoints), "\n")

		if (length(undo.breakpoints) < 1) {
			break
		}
		
		undo.breakpoints.shift <- undo.breakpoints + 1
		
		undo.df <- tempShort[undo.breakpoints, ]
		undo.df$seg.mean.diff <- abs(tempShort$seg.mean[undo.breakpoints] - tempShort$seg.mean[undo.breakpoints.shift])

		min.index <- which.min(undo.df$seg.mean.diff)
		
		leftIndex <- undo.df$segnum[min.index]
		rightIndex <- leftIndex + 1

		cat("sdundo.all left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"], tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"], tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
		cat("sdundo.all right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"], tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"], tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

		tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
		tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
		tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
		tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
		tempShort <- tempShort[-rightIndex, ]
		tempShort$segnum <- seq(1:nrow(tempShort))

	}
	
	return(tempShort)
	
}


cbs.segment01 <- function(indir, outdir, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
	#thisRatio$log2lowratio <- log(thisRatio$lowratio, base=2) 

	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), gc$chrom.arm, thisRatio$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	sortcol <- thisShort$chrom
	sortcol <- gsub("chr", "", sortcol)
	sortcol <- gsub("p", "", sortcol)
	sortcol <- gsub("q", "", sortcol)
	thisShort <- thisShort[order(as.numeric(sortcol)), ]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	cbs.long <- m[, 1]

	#####  NEW STUFF  also check min.width=2 above
	
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.short.cbs.txt", sep=""), quote=F, row.names=F) 

	workShort <- thisShort
	workShort$segnum <- 0
	workShort$seg.start <- 0
	workShort$seg.end <- 0
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		workShort$seg.start[i] <- thisStart
		workShort$seg.end[i] <- thisEnd
		workShort$segnum[i] <- i
		prevEnd = thisEnd
	}

	discardSegments <- TRUE
	while (discardSegments) {
		orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
		if (orderShort[1, "num.mark"] < min.width) {
			workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatio, undo.SD)
		} else {
			discardSegments <- FALSE
		}
	}

	workShort <- sdundo.all(workShort, thisRatio, undo.SD)
	thisShort <- workShort 
	
	#####  END NEW STUFF
	
	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	thisRatio$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatio$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatio$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatio$ratio.quantal <- thisRatio$lowratio * thisMultiplier
	thisRatio$seg.quantal <- thisRatio$seg.mean.LOWESS * thisMultiplier
	
	thisRatio$cbs.seg <- cbs.long
	thisRatio$cbs.seg.quantal <- cbs.long * thisMultiplier
	
	thisQuantalStats <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".50k.wg.png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".50k.wg.quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatio$abspos, y=thisRatio$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatio$abspos, y=thisRatio$ratio.quantal, col="#CCCCCC")
	points(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	lines(x=thisRatio$abspos, y=thisRatio$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStats, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.varbin.short.txt", sep=""), quote=F, row.names=F) 


	#bad <- read.table("/mnt/wigclust5/data/safe/kendall/badbins01/hg19.badbins.50k.txt", header=F, as.is=T, stringsAsFactors=F)
	bad <- read.table("/home/cnv/CNV_dataprocessing/badbins01/hg19.badbins.50k.txt", header=F, as.is=T, stringsAsFactors=F)

	thisRatioNobig <- thisRatio[-bad[, 1], ]

	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobig$lowratio, base=2), gc$chrom.arm[-bad[, 1]], thisRatioNobig$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	sortcol <- thisShort$chrom
	sortcol <- gsub("chr", "", sortcol)
	sortcol <- gsub("p", "", sortcol)
	sortcol <- gsub("q", "", sortcol)
	thisShort <- thisShort[order(as.numeric(sortcol)), ]

	m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	cbs.long.nobad <- m[, 1]

	#####  NEW STUFF  also check min.width=2 above

	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.short.cbs.txt", sep=""), quote=F, row.names=F) 

	workShort <- thisShort
	workShort$segnum <- 0
	workShort$seg.start <- 0
	workShort$seg.end <- 0
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		workShort$seg.start[i] <- thisStart
		workShort$seg.end[i] <- thisEnd
		workShort$segnum[i] <- i
		prevEnd = thisEnd
	}

	discardSegments <- TRUE
	while (discardSegments) {
		orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
		if (orderShort[1, "num.mark"] < min.width) {
			workShort <- remove.segment(workShort, orderShort[1, "segnum"], thisRatioNobig, undo.SD)
		} else {
			discardSegments <- FALSE
		}
	}

	workShort <- sdundo.all(workShort, thisRatioNobig, undo.SD)
	thisShort <- workShort

	#####  END NEW STUFF
	

	m <- matrix(data=0, nrow=nrow(thisRatioNobig), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatioNobig$seg.mean.LOWESS <- m[, 1]

	thisGrid <- seq(1.5, 5.5, by=0.05)
	thisOuter <- thisRatioNobig$seg.mean.LOWESS %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	thisShredded <- length(which(thisRatioNobig$seg.mean.LOWESS[which(chrom.numeric < 23)] < 0.1)) / length(which(chrom.numeric < 23))

	thisRatioNobig$ratio.quantal <- thisRatioNobig$lowratio * thisMultiplier
	thisRatioNobig$seg.quantal <- thisRatioNobig$seg.mean.LOWESS * thisMultiplier

	thisRatioNobig$cbs.seg <- cbs.long.nobad
	thisRatioNobig$cbs.seg.quantal <- cbs.long.nobad * thisMultiplier
	
	thisQuantalStatsNobad <- data.frame("ploidy"=thisMultiplier, "error"=thisError, "shredded"=thisShredded)
	
	chr <- thisRatioNobig$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	# Use the vlines abspos positions from above.  Because these start at
	# the third bin on the acrocentric chromosomes the vlines end up to
	# the right of the centromere rather than the left which is wrong.
	#vlines <- c(1, thisRatioNobig$abspos[which(chr != chr.shift) + 1], thisRatioNobig$abspos[nrow(thisRatioNobig)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

	png(paste(outdir, "/", sample.name, ".50k.wg.nobad.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$lowratio, col="#CCCCCC")
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA")
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.mean.LOWESS, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	hlines <- c(1, 2, 3, 4, 5, 6)

	png(paste(outdir, "/", sample.name, ".50k.wg.nobad.quantal.png", sep=""), height=800, width=1200)
	plot(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC")
	axis(1, at=x.at, labels=x.labels)
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$ratio.quantal, col="#CCCCCC")
	points(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	lines(x=thisRatioNobig$abspos, y=thisRatioNobig$seg.quantal, col="#0000AA")
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()

	write.table(thisQuantalStatsNobad, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.quantal.stats.txt", sep=""), quote=F, row.names=F) 
	write.table(thisRatioNobig, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.50k.k50.nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 


}

cbs.segment01(indir="/mnt/CNV/processed", outdir="/mnt/CNV/processed", varbin.gc="/home/cnv/CNV_dataprocessing/chrom.arm/varbin.gc.content.50k.bowtie.k50.hg19.txt", varbin.data="A832_4.varbin.50k.txt", sample.name="A832_4", alt.sample.name="20201028 lane x A832_4 A832 4 bcx x well x", alpha=0.02, nperm=1000, undo.SD=0.5, min.width=3)
