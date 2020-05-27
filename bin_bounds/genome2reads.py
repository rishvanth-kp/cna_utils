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
    for line in open(genomeFile):
        line = line.strip()
        if (line[0] == '>'):
            chrName = line[1:]
            genome.append([chrName, ''])
        else:
            genome[-1][1] += line
    return genome

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-genome', dest='genomeFile',
                        required=True, help='genome in fasta format')
    parser.add_argument('-l', '--read-len', dest='readLen',
                        type=int, required=True, help='read length')
    parser.add_argument('-n', '--reads-in-file', dest='readsInFile',
                        type=int, default=1000000, 
                        help='number of reads per file (default: 100M)')
    parser.add_argument('-o', '--outfile-prefix', dest='outfilePrefix',
                        required=True, help='output file name prefix')
    args = parser.parse_args()
    

    genome = readGenome(args.genomeFile)
    print(genome)

if __name__ == '__main__':
    main()
