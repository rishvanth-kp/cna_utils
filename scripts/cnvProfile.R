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
gc.smooth <- function(gc, ratio, lowess.f=0.05) {
  lowess.line <- lowess(gc, log2(ratio), f=lowess.f)
  lowess.points <- approx(lowess.line$x, lowess.line$y, gc)
  return(2^(log2(ratio) - lowess.points$y))
}


##
CBSsegment <- function(bin.counts, gc, min.width, seed, alpha,
                       n.perm, undo.sd, sample.name) {

  ## 4 colunm bed file. col 4: bin count
  bin.counts <- read.table(bin.counts)
  names(bin.counts) <- c("chr", "start", "end", "count")
  ## 5 column bed file. col 4: GC content, col5: chrom arm
  gc <- read.table(gc)
  names(gc) <- c("chr", "start", "end", "gc", "chr.arm")

  ## merge the bin counts and gc bed files. The result is order in
  ## bin.counts is preserved
  seg <- merge(bin.counts, gc, by=c("chr", "start", "end"), sort=FALSE)

  seg$ratio <- (seg$count + 1) / mean(seg$count + 1)
  seg$lowess.ratio <- gc.smooth(seg$gc, seg$ratio)
  print(head(seg))

  set.seed(seed)
  cbs.seg <- smooth.CNA(CNA(log2(seg$lowess.ratio), seg$chr.arm,
                          seg$start, data.type="logratio",
                          sampleid=sample.name))
  cbs.seg <- segment(cbs.seg, alpha=alpha, nperm=n.perm,
                undo.splits="sdundo", undo.SD=undo.sd, min.width=2)
  print(cbs.seg[[2]])

  return(seg)
}

main <- function() {

  parser <- OptionParser()
  parser <- add_option(parser, c("-b", "--bincounts"),
              help="bin counts bed file [Required]")
  parser <- add_option(parser, c("-g", "--gc"),
              help="GC content bed file [Required]")
  parser <- add_option(parser, c("-n", "--samplename"),
              help="Sample name [Required]")
  parser <- add_option(parser, c("-w", "--minwidth"), default=3,
              help="CBS min segment width [Default: %default]")
  parser <- add_option(parser, c("-s", "--seed"), default=25,
              help="RNG seed for CBS [Default: %default]")
  parser <- add_option(parser, c("-a", "--alpha"), default=0.02,
              help="CBS segment alpha [Default: %default]")
  parser <- add_option(parser, c("-p", "--nperm"), default=1000,
              help="CBS number of permutations [Default: %default]")
  parser <- add_option(parser, c("-u", "--undosd"), default=0.5,
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
                min.width=opt$minwidth, seed=opt$seed,
                alpha=opt$alpha, n.perm=opt$nperm, undo.sd=opt$undosd,
                sample.name=opt$samplename)

  write.table(cbs.seg, "test.txt", quote=FALSE)
}


main()
