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

n.gt <- function(x, gt.val) {
  return (sum(x > gt.val))
}

main <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("determineBadBins.R <list-of-segs.txt> <out-bad-bins.txt>,
          call.=FALSE")
  }

  sample.files <- read.table(args[1], header = F)
  names(sample.files) <- c("files")
  outfile <- args[2]

  first <- TRUE
  ratios <- data.frame()
  n.samples <- 0
  for (i in sample.files$files) {
    ## first 3 colums contain chr, start, and end
    ## column 8 has the lowratio
    sample <- read.table(i, header = TRUE)
    sample.name <- sprintf("sample.%d", n.samples)
    if (first) {
      ratios <- sample[, c(1, 2, 3, 8)]
      names(ratios) <- c("chrom", "start", "end", sample.name)
      first <- FALSE
    }
    else {
      ratios$lowratio <- sample[, 8]
      names(ratios)[ncol(ratios)] <- sample.name
    }
    n.samples <- n.samples + 1
  }

  n.bins <- nrow(ratios)

  ## exclude sex chromosomes
  ratios <- ratios[!ratios$chrom == 23,]
  ratios <- ratios[!ratios$chrom == 24,]

  ## standardize the bin ratios using the mean and sd of each
  ## chromosome.
  for (i in seq(1, 22)) {
    ratios[ratios$chrom == i, 4:ncol(ratios)] <-
      (ratios[ratios$chrom == i, 4:ncol(ratios)] -
      mean(as.matrix(ratios[ratios$chrom == i, 4:ncol(ratios)]))) /
      sd(as.matrix(ratios[ratios$chrom == i, 4:ncol(ratios)]))
  }

  ## cutoff of calling bad-bins is based on the p-value corresponding
  ## to 1/(number of bins)
  cut.off <- abs(qnorm(1 / n.bins))
  # cut.off <- 2

  ## Designate as a bad bin if the median of the ratio over all the samples
  ## is above the cutoff
  # ratios$median <- apply(as.matrix(ratios[, 4:ncol(ratios)]), 1,  median)
  # bad.bins <- ratios[ratios$median >= cut.off, 1:3]
  ratios$median <- apply(as.matrix(ratios[, 4:ncol(ratios)]), 1, n.gt,
                          gt.val = cut.off)
  bad.bins <- ratios[ratios$median >= (n.samples/2), 1:3]


  write.table(bad.bins, outfile, sep = "\t", quote = FALSE,
              col.names = FALSE)
}

main()
