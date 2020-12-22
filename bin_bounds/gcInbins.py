#!/usr/bin/env python3

# Copyright (C) 2019 Rishvanth Prabakar
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

def readGenome(genomeFile):
    genome = {}
    first = True
    chrName = ''
    chrSeq = ''
    for line in open(genomeFile):
        line = line.strip()
        if (line[0] == '>'):
            if (not first):
                genome[chrName] = chrSeq
            chrName = line[1:]
            chrSeq = ''
            first = False
        else:
            chrSeq += line.upper()
    genome[chrName] = chrSeq

    return genome


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', dest = 'genomeFile',
                        required = True, help = 'genome in fasta format')
    parser.add_argument('-b', '--bin-bounds', dest = 'binBounds',
                        required = True, help = 'bin boundaries BED file')
    parser.add_argument('-c', '--chrom-arm', dest = 'chromArms',
                        required = True, help = 'start coordinates of q-arm')
    parser.add_argument('-o', '--outfile', dest = 'outfile',
                        required = True, help = 'output file name')
    args = parser.parse_args()


    genome = readGenome(args.genomeFile)

    chromArm = {}
    chromArmFile = open(args.chromArms, 'r')
    chromArmFile.readline()
    for line in chromArmFile:
        line = line.strip().split()
        chromArm[line[0]] = int(line[2])

    out = open(args.outfile, 'w')
    for line in open(args.binBounds):
        line = line.strip().split()
        ch = line[0]
        start = int(line[1])
        end = int(line[2])
        gCount = genome[ch].count('G', start, end)
        cCount = genome[ch].count('C', start, end)
        gcContent = (gCount + cCount) / (end - start)
        if (start < chromArm[ch]):
            arm = ch + 'p'
        else:
            arm = ch + 'q'
        print('%s\t%d\t%d\t%f\t%s' % 
                (ch, start, end, gcContent, arm), 
                file = out)

    out.close()


if __name__ == '__main__':
    main()
