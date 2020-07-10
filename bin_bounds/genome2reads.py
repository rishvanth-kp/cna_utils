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
    genome = []
    genomeLen = 0
    for line in open(genomeFile):
        line = line.strip()
        if (line[0] == '>'):
            chrName = line[1:]
            genome.append([chrName, ''])
        else:
            genome[-1][1] += line
            genomeLen += len(line)
    return genome, genomeLen


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-genome', dest='genomeFile',
                        required=True, help='genome in fasta format')
    parser.add_argument('-l', '--read-len', dest='readLen',
                        type=int, required=True, help='read length')
    parser.add_argument('-n', '--reads-in-file', dest='readsInFile',
                        type=int, default=25000000, required=False,
                        help='number of reads per file (default: 25M)')
    parser.add_argument('-o', '--outfile-prefix', dest='outfilePrefix',
                        required=True, help='output file name prefix')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help='print progress to stderr')
    args = parser.parse_args()

    if (args.VERBOSE):
        print('[READING GENOME]', file = sys.stderr)
    genome, genomeLen = readGenome(args.genomeFile)
    if (args.VERBOSE):
        print('\tChromosomes: %d' % (len(genome)), file=sys.stderr)
        print('\tGenome len: %d' % (genomeLen), file=sys.stderr)

    if (args.VERBOSE):
        print('[WRITING READS]', file = sys.stderr)
    fileCount = 0
    for chrs in genome:
        if (args.VERBOSE):
            print('\t%s: %d' % (chrs[0], len(chrs[1])))
        for i in range(0, len(chrs[1]) - args.readLen + 1):
            if (i % args.readsInFile == 0):
                fileName = '%s_%d.fa' % (args.outfilePrefix, fileCount)
                out = open(fileName, 'w')
                fileCount += 1
            readName = '%s:%d' % (chrs[0], i)
            readSeq = chrs[1][i:i+args.readLen]
            print('%s\n%s' % (readName, readSeq), file = out)


if __name__ == '__main__':
    main()
