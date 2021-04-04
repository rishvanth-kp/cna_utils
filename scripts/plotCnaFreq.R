#!/usr/bin/env Rscript

# Copyright (C) 2020 Rish Prabakar
#
# Authors: Rish Prabakar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

main <- function() {

  args = commandArgs(trailingOnly=TRUE)
  if (length(args) != 5) {
    stop ("plotCnaFreq.R <cna-csv> <ch-bin-start>
            <gain-threshold> <loss-threshold> <sample-name>", call.=FALSE)
  }

  ## csv file downloaded from the CNV app. First three columns have chrom,
  ## chrompos, and abosolute chompos. Followed by one column per sample 
  cna <- read.csv(args[1], na.strings="", header=T)
  # 3 column file with chr, bin.start, and bin.end (starting bin for chr)
  chr.bounds <- read.table(args[2], header=T)
  gain.thresh <- as.numeric(args[3])
  loss.thresh <- as.numeric(args[4])
  sample.name <- args[5]

  ## get rid of the first three columns
  cna <- t(cna[-1:-3])
  seg.col.start <- 1
  seg.col.end <- ncol(cna)
  print(dim(cna))

  n.samples <- length(cna[, 1])
  print(n.samples)

  amp <- c()
  for (i in seq(seg.col.start, seg.col.end)) {
    amp <- c(amp, sum((cna[, i]) >= gain.thresh))
  }
  amp <- amp / n.samples

  del <- c()
  for (i in seq(seg.col.start, seg.col.end)) {
    del <- c(del, sum((cna[, i]) <= loss.thresh))
  }
  del <- del / n.samples

  y.max <- max(amp) + 0.05
  y.min <- min(-del) - 0.05

  x.at <- (chr.bounds$bin.start + chr.bounds$bin.end) / 2
  x.labels <- chr.bounds$chr
  if (length(x.labels) == 24) {
    x.labels[23:24] <- c("X", "Y") 
  }
  y.at <- c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)
  y.labels <- c(1, 0.8, 0.6, 0.4, 0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)


  pdf(sprintf("%s.pdf", sample.name), width=6.5, height=3.5, useDingbats=F)
  par(pin=c(5,2.25))
  plot(amp, ylim=c(y.min, y.max), type="n", xaxs="i", yaxs="i", axes=FALSE,
        ann=FALSE)
  rect(chr.bounds$bin.start[seq(1, length(chr.bounds$bin.start), 2)], y.min,
       chr.bounds$bin.end[seq(1, length(chr.bounds$bin.end), 2)], y.max,
       col="lightgrey", border="NA")
  lines(amp, col="red")
  lines(-del, col="blue")
  abline(0, 0)
  axis(side=2, at=y.at, labels=y.labels)
  axis(side=1, at=x.at, labels=FALSE)
  mtext(text=x.labels, side=1, at=x.at, line=c(0.25,0.75), cex=0.7)
  mtext(text="Chromosome", side=1, line=2, cex=1.1)
  mtext(text="Frequency", side=2, line=2, cex=1.1)
  mtext(text=sprintf("%s (n=%d)", sample.name, n.samples), side=3, 
        line=1, cex=1.2)
  box()
  dev.off()

}

main ()
