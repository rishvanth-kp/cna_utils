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
    parser.add_argument('-i', '--input-bam', dest='bamFile',
                        required=True, help='mapped bam file')
    parser.add_argument('-q', '--min-mapq', dest='minMapq', type=int,
                        default=0, required=False, help='min MAPQ to consider')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='out file in BED format')
    args = parser.parse_args()


    chrStart = 0
    chrEnd = 0
    chrName = ''
    first = True;
    alignedReads = pysam.AlignmentFile(args.bamFile, 'rb')
    out = open(args.outFile, 'w')
    for aln in alignedReads:
        if (aln.mapping_quality >= args.minMapq and
            int(aln.query_name.split(':')[1]) == aln.reference_start):
            if first:
                chrName = aln.reference_name
                chrStart = aln.reference_start
                first = False
            chrEnd = aln.reference_start + 1
        else:
            if (chrStart != chrEnd):
                print('%s\t%d\t%d' % (chrName, chrStart, chrEnd),
                        file = out)
            chrStart = chrEnd
            first = True
    if (chrStart != chrEnd):
        print('%s\t%d\t%d' % (chrName, chrStart, chrEnd),
                file = out)

    alignedReads.close()
    out.close()

if __name__ == '__main__':
    main()
