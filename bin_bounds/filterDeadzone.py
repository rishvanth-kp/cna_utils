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
import bisect
import pysam


def isDead(chrom, pos, deadzones):
    if chrom in deadzones:
        offset = bisect.bisect_left(deadzones[chrom], (pos, sys.maxsize))
        if (offset > 0):
            if (pos >= deadzones[chrom][offset-1][0] and
                pos < deadzones[chrom][offset-1][1]):
                return True
    return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-bam', dest='bamFile',
                        required=True, help='mapped bam file')
    parser.add_argument('-b', '--dz-bed', dest='bedFile',
                        required=True, help='deadzone bed file')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='out file in SAM/BAM format')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help='print progress to stderr')
    args = parser.parse_args()


    if (args.VERBOSE):
        print('[READING DEADZONES]', file=sys.stderr)
    deadzones = {}
    for line in open(args.bedFile):
        chrom, start, end = line.strip().split()
        if chrom in deadzones:
            deadzones[chrom].append((int(start), int(end)))
        else:
            deadzones[chrom] = [(int(start), int(end))]
    for vals in deadzones.values():
        vals.sort()

    if (args.VERBOSE):
        print('[FILTERING ALIGNMENTS]', file=sys.stderr)
    totalMaps = 0
    passedMaps = 0
    writeFlags = 'w'
    if (args.outFile[-3:] == 'bam'):
        writeFlags += 'b'
    aligned = pysam.AlignmentFile(args.bamFile, 'r')
    filtered = pysam.AlignmentFile(args.outFile, writeFlags, template=aligned)
    for aln in aligned:
        alnChrom = aln.reference_name
        alnPos = aln.reference_start
        totalMaps += 1
        if (not isDead(alnChrom, alnPos, deadzones)):
            filtered.write(aln)
            passedMaps += 1
    aligned.close()
    filtered.close()

    if (args.VERBOSE):
        filteredMaps = totalMaps - passedMaps
        passedFrac = (passedMaps / totalMaps) * 100
        filteredFrac = (filteredMaps / totalMaps) * 100
        print('[FILTERED STATS]')
        print('\tTotal alignments: %d' % (totalMaps))
        print('\tPassed alignments: %d (%f%%)' %
                (passedMaps, passedFrac))
        print('\tFiltered alignments: %d (%f%%)' %
                (filteredMaps, filteredFrac))

if __name__ == '__main__':
   main()
