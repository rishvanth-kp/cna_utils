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
  if (length(args) != 3) {
    stop("binBoundsAfterExclude.R <bin-bounds.bed> <exclude-bins.txt>
                      <outfile.bed>", call.=FALSE)
  }

  bin.bounds <- read.table(args[1], header=FALSE)
  names(bin.bounds) <- c("chr", "start", "end")
  exclude <- read.table(args[2], header=FALSE)
  names(exclude) <- c("bins")
  outfile <- args[3]

  # collapse the bin boundaries
  i <- 1
  exclude.collapsed <- c()
  while (i <= nrow(exclude)) {
    j <- 1
    while(i != nrow(exclude) & 
            ((exclude$bins[i + j] - exclude$bins[i]) == j)) {
      j <- j + 1
    }
    exclude.collapsed <- rbind(exclude.collapsed, 
                                c(exclude$bins[i], j))
    i <- i + j
  }
 
  exclude.collapsed <- as.data.frame(exclude.collapsed) 
  names(exclude.collapsed) <- c("bin", "length")


  # Fix the boundaries of bins that are adjacent to the excluded bins
  # If the excluded bin occures in the middle of the chromosone, then
  #   extend the end of the previous bin to the end of the excluded bin.
  # If the excluded bin in the first one on the chromosome, then
  #   set the start of the next bin to start of the excluded bin (which
  #   should always be 0).
  
  for (i in 1:nrow(exclude.collapsed)) {
    rm.bin <- exclude.collapsed[i,]$bin
    rm.len <- exclude.collapsed[i,]$len
    if (bin.bounds[rm.bin,]$start == 0) {
      bin.bounds[rm.bin+rm.len,]$start = bin.bounds[rm.bin,]$start
    }
    else {
      bin.bounds[rm.bin-1,]$end = bin.bounds[rm.bin+rm.len-1,]$end
    }
  }


  # remove the excluded bins
  bin.bounds.excluded <- bin.bounds[-exclude$bins,]
  print(dim(bin.bounds.excluded))

  write.table(bin.bounds.excluded, outfile, sep="\t",
              row.names=FALSE, quote=FALSE)

}

main()
