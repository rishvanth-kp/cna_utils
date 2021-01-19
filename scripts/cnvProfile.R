#!/usr/bin/env Rscript

# cnvProfile.R: Adapted from "SRR054616.cbs.r" supplemetary program
# 12 of Baslan, Timour, et al. "Genome-wide copy number analysis of
# single cells." Nature protocols 7.6 (2012): 1024.
#
# Functions RemoveSegment and SDUndoAll are adapted from the
# modification to supplemetary program 12 by Jude Kendall.
#
# Copyright (C) 2020 CSI-Cancer, University of Southern California
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

library("optparse")
library("DNAcopy")

##
GCsmooth <- function(gc, ratio, lowess.f=0.05) {
  lowess.line <- lowess(gc, log2(ratio), f=lowess.f)
  lowess.points <- approx(lowess.line$x, lowess.line$y, gc)
  return(2^(log2(ratio) - lowess.points$y))
}


## The output of the segment function is ordered alphabatically by
## the chromosome name. This function converts it to chromosome order.
CBSsort <- function(seg) {
  chrom.numeric <- gsub("chr|p|q", "", seg$chrom)
  chrom.numeric <- gsub("X", 23, chrom.numeric)
  chrom.numeric <- gsub("Y", 24, chrom.numeric)
  chrom.numeric <- as.numeric(chrom.numeric)
  seg <- seg[order(chrom.numeric),]
  return(seg)
}


##
Short2LongSegments <- function(seg) {

  # print("Short2LongSegments")
  # print(sum(seg$num.mark))
  long.seg <- c()
  for (i in 1:nrow(seg)) {
    long.seg <- c(long.seg,
                    rep(seg$seg.mean[i], seg$num.mark[i]))
  }
  long.seg <- 2^long.seg
  return(long.seg)
}


##
RemoveSegment <- function(seg, bin.ratio, undo.sd, index) {
  print("REMOVING SEGMENT")
  print(index)
  print(seg[index,])

  append.left <- TRUE
  check.sd.undo <- FALSE
  if (index == 1) {
    append.left <- FALSE
  }
  else if (index == nrow(seg)) {
    append.left <- TRUE
  }
  else if (seg$chrom[index] != seg$chrom[index - 1]) {
    append.left <- FALSE
  }
  else if (seg$chrom[index] != seg$chrom[index + 1]) {
    append.left <- TRUE
  }
  else {
    check.sd.undo <- TRUE
    if (abs(seg$seg.mean[index - 1] - seg$seg.mean[index]) <
        abs(seg$seg.mean[index + 1] - seg$seg.mean[index])) {
      append.left <- TRUE
    }
    else {
      append.left <- FALSE
    }
  }

  append.index <- index + 1
  if (append.left) {
    append.index <- index - 1
  }

  print(append.left)
  print(check.sd.undo)
  print(append.index)

  if (append.left) {
    seg$loc.end[append.index] <- seg$loc.end[index]
    seg$end[append.index] <- seg$end[index]
  }
  else {
    seg$loc.start[append.index] <- seg$loc.start[index]
    seg$start[append.index] <- seg$start[index]
  }

  seg$num.mark[append.index] <- (seg$num.mark[append.index] +
                                  seg$num.mark[index])
  ## This is bad design. It makes the assumption that the chromosome order
  ## in seg and bin ratio are the same
  seg$seg.mean[append.index] <- mean(log2(bin.ratio[seg$start[append.index]:
                                          seg$end[append.index]]))
  seg <- seg[-index,]
  seg$num <- seq(1:nrow(seg))

  if (check.sd.undo) {
    left.index <- index - 1
    right.index <- index

    print(seg[left.index,])
    print(seg[right.index,])

    bin.ratio.sd <- mad(diff(bin.ratio)) / sqrt(2)
    if (abs(seg$seg.mean[left.index] - seg$seg.mean[right.index]) <
          (bin.ratio.sd * undo.sd)) {
      print("UNDO SD IN REMOVE SEGMENT")
      seg$loc.end[left.index] <- seg$loc.end[right.index]
      seg$end[left.index] <- seg$end[right.index]
      print("FUCK")
      print(seg[left.index,])
      seg$num.mark[left.index] <- seg$num.mark[left.index] +
                                    seg$num.mark[right.index]
      seg$seg.mean[left.index] <- mean(log2(bin.ratio[seg$start[left.index]:
                                            seg$end[right.index]]))
      seg <- seg[-right.index,]
      seg$num <- seq(1:nrow(seg))
    }
  }

  print(seg)
  return(seg)
}


##
SegmentsUndoSD <- function(seg, bin.ratio, undo.sd) {

  bin.ratio.sd <- mad(diff(bin.ratio)) / sqrt(2)

  chrom <- seg$chrom
  chrom.shift <- c(seg$chrom[-1], seg$chrom[1])
  breakpoints <- which(chrom == chrom.shift)

  undo.breakpoints <- breakpoints[which(abs(seg$seg.mean[breakpoints] -
                                  seg$seg.mean[breakpoints + 1]) <
                                  bin.ratio.sd * undo.sd)]
  # print("REMOVING SEGMENT IN UNDO SD ALL")
  # print(breakpoints)
  # print(undo.breakpoints)

  while(length(undo.breakpoints) >= 1) {
    undo.df <- seg[undo.breakpoints,]
    undo.df$seg.mean.diff <- abs(seg$seg.mean[undo.breakpoints] -
                                 seg$seg.mean[undo.breakpoints + 1])

    # print(undo.df)
    left.index <- undo.df$num[which.min(undo.df$seg.mean.diff)]
    right.index <- left.index + 1
    # print(left.index)
    seg$loc.end[left.index] <- seg$loc.end[right.index]
    seg$end[left.index] <- seg$end[right.index]
    seg$num.mark[left.index] <- seg$num.mark[left.index] +
                                  seg$num.mark[right.index]
    seg$seg.mean[left.index] <- mean(log2(bin.ratio[seg$start[left.index]:
                                          seg$end[right.index]]))
    seg <- seg[-right.index,]
    seg$num <- seq(1:nrow(seg))

    chrom <- seg$chrom
    chrom.shift <- c(seg$chrom[-1], seg$chrom[1])
    breakpoints <- which(chrom == chrom.shift)
    undo.breakpoints <- breakpoints[which(abs(seg$seg.mean[breakpoints] -
                                    seg$seg.mean[breakpoints + 1]) <
                                    bin.ratio.sd * undo.sd)]
    # print(breakpoints)
    # print(undo.breakpoints)
  }

  return(seg)
}


##
RemoveShortSegments <- function(seg, bin.ratio , min.width, undo.sd) {

  seg$end <- cumsum(seg$num.mark)
  seg$start <- seg$end - seg$num.mark + 1
  seg$num <- seq(1:nrow(seg))

  seg <- seg[,c("ID", "chrom", "loc.start", "loc.end", "num.mark",
                "start", "end", "num", "seg.mean")]

  print(seg)

  while (min(seg$num.mark) < min.width) {
    seg <- RemoveSegment(seg, bin.ratio, undo.sd,
              seg$num[order(seg$num.mark, abs(seg$seg.mean))[1]])
  }

  seg <- SegmentsUndoSD(seg, bin.ratio, undo.sd)

  return(seg)
}


##
MergeAcrocentric <- function(seg, min.width) {

  for (i in unique(seg$chr)) {
    arms <- unique(seg[seg$chr == i,]$chr.arm)
    if (length(arms) == 2) {
      if (nrow(seg[seg$chr.arm == arms[1],]) < min.width |
          nrow(seg[seg$chr.arm == arms[2],]) < min.width) {
        seg[seg$chr == i,]$chr.arm = i
      }
    }
  }
  return(seg)
}


##
CBSquantal <- function(seg, min.ploidy, max.ploidy) {

  ploidy <- seq(min.ploidy, max.ploidy, by=0.05)
  seg.cn <- outer(seg$seg.mean, ploidy)
  cn.error <- (seg.cn - round(seg.cn)) ^ 2
  cn.error <- colSums(cn.error)
  ploidy <- ploidy[which.min(cn.error)]
  cn.error <- min(cn.error)

  seg$lowess.ratio.quantal <- seg$lowess.ratio * ploidy
  seg$seg.mean.quantal <- seg$seg.mean * ploidy
  # print(ploidy)

  return(seg)
}

##
CBSsegment <- function(bin.counts, gc, bad.bins=NULL, min.width, min.ploidy,
                       max.ploidy, seed, alpha, n.perm, undo.sd, sample.name) {

  ## 4 colunm bed file. col 4: bin count
  bin.counts <- read.table(bin.counts)
  names(bin.counts) <- c("chr", "start", "end", "count")
  ## 5 column bed file. col 4: GC content, col5: chrom arm
  gc <- read.table(gc)
  names(gc) <- c("chr", "start", "end", "gc", "chr.arm")

  ## merge the bin counts and gc bed files. The order in
  ## bin.counts is preserved
  seg <- merge(bin.counts, gc, by=c("chr", "start", "end"), sort=FALSE)

  ## Merge p and q chromosome arms if a chromosome arm is shorter than
  ## the minimum segment width
  seg <- MergeAcrocentric(seg, min.width)

  seg$ratio <- (seg$count + 1) / mean(seg$count + 1)
  seg$lowess.ratio <- GCsmooth(seg$gc, seg$ratio)

  if (!is.null(bad.bins)) {
    # print("Removing bad bins")
    bad.bins <- read.table(bad.bins)
    names(bad.bins) <- c("chr", "start", "end")
    # print(bad.bins)
    # print(head(seg[,1:3]))
    bin.match <- match(paste(seg[,1], seg[,2], sep="_"),
                       paste(bad.bins[,1], bad.bins[,2], sep="_"),
                       nomatch=FALSE)
    seg <- seg[!bin.match,]
  }

  set.seed(seed)
  seg.short <- smooth.CNA(CNA(log2(seg$lowess.ratio), seg$chr.arm,
                          seg$start, data.type="logratio",
                          sampleid=sample.name))
  seg.short <- segment(seg.short, alpha=alpha, nperm=n.perm,
                undo.splits="sdundo", undo.SD=undo.sd, min.width=2)
  seg.short <- seg.short[[2]]

  seg.short <- CBSsort(seg.short)
  seg.short <- RemoveShortSegments(seg.short, seg$lowess.ratio,
                min.width, undo.sd)

  seg$seg.mean <- Short2LongSegments(seg.short)
  seg <- CBSquantal(seg, min.ploidy, max.ploidy)

  return(list(long.seg=seg, short.seg=seg.short))
}

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-b", "--bincounts"),
              help="bin counts bed file [Required]")
  parser <- add_option(parser, c("-g", "--gc"),
              help="GC content bed file [Required]")
  parser <- add_option(parser, c("-e", "--badbins"),
              help="Bad bins to exclude bed file [Optional]")
  parser <- add_option(parser, c("-n", "--samplename"),
              help="Sample name [Required]")
  parser <- add_option(parser, c("-w", "--minwidth"), default=3,
              help="CBS min segment width [Default: %default]")
  parser <- add_option(parser, c("--minploidy"), default=1.5,
              help="Min. ploidy for CN estimation [Default: %default]")
  parser <- add_option(parser, c("--maxploidy"), default=5.5,
              help="Max. ploidy for CN estimation [Default: %default]")
  parser <- add_option(parser, c("--seed"), default=25,
              help="RNG seed for CBS [Default: %default]")
  parser <- add_option(parser, c("--alpha"), default=0.02,
              help="CBS segment alpha [Default: %default]")
  parser <- add_option(parser, c("--nperm"), default=1000,
              help="CBS number of permutations [Default: %default]")
  parser <- add_option(parser, c("--undosd"), default=0.5,
              help="CBS undo SD [Default: %default]")
  parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
              default=FALSE, help="Print extra output [FALSE]")
  opt <- parse_args(parser)

  if (is.null(opt$bincounts) | is.null(opt$gc) |
      is.null(opt$samplename)) {
    print_help(parser)
    quit(status=1)
  }

  cbs.seg <- CBSsegment(bin.counts=opt$bincounts, gc=opt$gc,
                min.width=opt$minwidth, min.ploidy=opt$minploidy,
                max.ploidy=opt$maxploidy, seed=opt$seed,
                alpha=opt$alpha, n.perm=opt$nperm, undo.sd=opt$undosd,
                sample.name=opt$samplename)

  write.table(cbs.seg$short.seg, sprintf("%s_short_seg.txt", opt$samplename),
    row.names=FALSE, sep="\t", quote=FALSE)
  write.table(cbs.seg$long.seg, sprintf("%s_seg.txt", opt$samplename),
    row.names=FALSE, sep="\t", quote=FALSE)

  if (!is.null(opt$badbins)) {
    cbs.seg.nobad <- CBSsegment(bin.counts=opt$bincounts, gc=opt$gc,
                        bad.bins=opt$badbins, min.width=opt$minwidth,
                        min.ploidy=opt$minploidy, max.ploidy=opt$maxploidy,
                        seed=opt$seed, alpha=opt$alpha, n.perm=opt$nperm,
                        undo.sd=opt$undosd, sample.name=opt$samplename)

    write.table(cbs.seg.nobad$short.seg, sprintf("%s_short_seg_nobad.txt",
      opt$samplename), row.names=FALSE, sep="\t", quote=FALSE)
    write.table(cbs.seg.nobad$long.seg, sprintf("%s_seg_nobad.txt",
      opt$samplename), row.names=FALSE, sep="\t", quote=FALSE)
  }

}

main()
