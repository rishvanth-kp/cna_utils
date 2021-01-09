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
Short2LongSegments <- function() {

}

##
RemoveSegment <- function(seg, bin.ratio, undo.sd, index) {
  print("Removing segment")
  print(index)
  print(seg[index,])
  print(dim(seg))
 
  append.left <- TRUE
  check.sd.undo <- FALSE
  ## TODO: HOLY FUCK!!! NEED TO GET RID OF THE CODE BELOW ASAP
  if (index == 1) {
    append.left <- FALSE
  }
  else {
    if (index == nrow(seg)) {
      append.left <- TRUE
    }
    else {
      right.index <- index + 1
      left.index <- index - 1
    
      if (seg$chrom[right.index] != seg$chrom[index]) {
        append.left <- TRUE
      }
      else {
        if (seg$chrom[left.index] != seg$chrom[index]) {
          append.left <- FALSE
        }
        else {
          if (abs(seg$seg.mean[left.index] - seg$seg.mean[index]) < 
              abs(seg$seg.mean[right.index] - seg$seg.mean[index])) {
            append.left <- TRUE
            check.sd.undo <- TRUE
          }
          else {
            append.left <- FALSE
            check.sd.undo <- TRUE   
          }
        }
      }
    }
  }


  append.index <- index + 1
  if (append.left) {
    append.index <- index - 1
  } 

  print(append.index)
  print(append.left)
 
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
  ## TODO: Add code to verify that that bin.ratio and seg are in the same 
  ## chromosome order
  seg$seg.mean[append.index] <- mean(log2(bin.ratio[seg$start[append.index]:
                                          seg$end[append.index]]))
  seg <- seg[-index,]
  seg$num <- seq(1:nrow(seg))

  if (check.sd.undo) {
    left.index <- index - 1
    right.index <- index

    bin.ratio.sd <- mad(diff(bin.ratio)) / sqrt(2)
    if (abs(seg$seg.mean[left.index] - seg$seg.mean[right.index]) < 
          (bin.ratio.sd * undo.sd)) {
      print("UNDO SD IN REMOVE SEGMENT") 
      seg$loc.end[left.index] <- seg$loc.end[right.index]
      seg$seg.end[left.index] <- seg$seg.end[right.index]
      seg$num.mark[left.index] <- seg$num.mark[left.index] + 
                                    seg$num.mark[right.index]
      seg$seg.mean[left.index] <- mean(log2(bin.ratio[seg$start[left.index]:
                                            seg$end[right.index]]))
      seg <- seg[-right.index,]
      seg$num <- seq(1:nrow(seg))
    }
  }

  return(seg)
}

SegmentsUndoSD <- function(seg, bin.ratio, undo.sd) {
  
  bin.ratio.sd <- mad(diff(bin.ratio)) / sqrt(2)

  chrom <- seg$chrom
  chrom.shift <- c(seg$chrom[-1], seg$chrom[1])
  breakpoints <- which(chrom == chrom.shift)

  undo.breakpoints <- breakpoints[which(abs(seg$seg.mean[breakpoints] - 
                                  seg$seg.mean[breakpoints + 1]) <
                                  bin.ratio.sd * undo.sd)]
  print(breakpoints) 
  print(undo.breakpoints)
  
  while(length(undo.breakpoints) >= 1) {
    print("REMOVING SEGMENT IN UNDO SD ALL")
    undo.df <- seg[undo.breakpoints,]
    undo.df$seg.mean.diff <- abs(seg$seg.mean[undo.breakpoints] - 
                                 seg$seg.mean[undo.breakpoints + 1])

    print(undo.df)
    left.index <- undo.df$num[which.min(undo.df$seg.mean.diff)]
    right.index <- left.index + 1
    print(left.index)
    seg$loc.end[left.index] <- seg$loc.end[right.index]
    seg$seg.end[left.index] <- seg$seg.end[right.index]
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
    print(breakpoints) 
    print(undo.breakpoints)
  }
  
     

  return(seg)
}

## 
RemoveShortSegments <- function(seg, bin.ratio , min.width, undo.sd) {
 
  ## TODO: The code below can be replaced with 'cumsum' and 'seq' 
  seg$num <- 0
  seg$start <- 0
  seg$end <- 0
  prev.end <- 0
  for (i in 1:nrow(seg)) {
    start <- prev.end + 1
    end <- prev.end + seg$num.mark[i]
    seg$start[i] <- start
    seg$end[i] <- end
    seg$num[i] <- i 
    prev.end <- end
  }
 
  while (min(seg$num.mark) < min.width) {
    seg <- RemoveSegment(seg, bin.ratio, undo.sd, 
              seg$num[order(seg$num.mark, abs(seg$seg.mean))[1]])
  }
 
  seg <- SegmentsUndoSD(seg, bin.ratio, undo.sd)

  return(seg)
} 


##
MergeAcrocentric <- function(seg) {

  return(seg)
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

  ## merge the bin counts and gc bed files. The order in
  ## bin.counts is preserved
  seg <- merge(bin.counts, gc, by=c("chr", "start", "end"), sort=FALSE)

  seg$ratio <- (seg$count + 1) / mean(seg$count + 1)
  seg$lowess.ratio <- GCsmooth(seg$gc, seg$ratio)

  set.seed(seed)
  cbs.seg <- smooth.CNA(CNA(log2(seg$lowess.ratio), seg$chr.arm,
                          seg$start, data.type="logratio",
                          sampleid=sample.name))
  cbs.seg <- segment(cbs.seg, alpha=alpha, nperm=n.perm,
                undo.splits="sdundo", undo.SD=undo.sd, min.width=2)
  cbs.seg <- cbs.seg[[2]]

  ## TODO: need to deal with acrocentric chromosomes. At present, they 
  ## generate bins of width 1. And RemoveShortSegments sometimes merges 
  ## these with the wrong chomosomes.

  cbs.seg <- CBSsort(cbs.seg)
  print(cbs.seg)
  cbs.seg <- RemoveShortSegments(cbs.seg, seg$lowess.ratio, 
                min.width, undo.sd)
  print(cbs.seg)

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

  write.table(cbs.seg, sprintf("%s_seg.txt", opt$samplename), 
    quote=FALSE)
}


main()
