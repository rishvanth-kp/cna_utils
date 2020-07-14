#!/usr/bin/env python3

# Copyright (C) 2020 Rishvanth Prabakar
#
# Authors: Rish Prabakar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import sys, argparse
import pysam

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-sam', dest='samFile',
                        required=True, help='mapped sam file')
    parser.add_argument('-q', '--min-mapq', dest='minMapq', type=int, 
                        default=0, required=False, help='min MAPQ to consider')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='out file')
    args = parser.parse_args()

    alignmentLen = {}
    totalMapped = 0  
    mappedReads = pysam.AlignmentFile(args.samFile, 'r')
    for read in mappedReads:
        if read.mapping_quality >= args.minMapq:
            totalMapped += 1
            alnLen = read.query_alignment_length
            if alnLen in alignmentLen:
                alignmentLen[alnLen] += 1
            else:
                alignmentLen[alnLen] = 1
    mappedReads.close()

    out = open(args.outFile, 'w')
    for length, count in sorted(alignmentLen.items()):
        print('%d\t%d\t%f' % 
                (length, 
                count,
                count/totalMapped), 
                file=out)
    out.close()

if __name__ == '__main__':
    main()
