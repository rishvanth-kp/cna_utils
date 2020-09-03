#!/usr/bin/env python3

# Copyright (C) 2020 Kuhn-Hicks Lab, University of Southern California
#
# Authors: Rish Prabakar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, argparse
import gzip
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-fastq', dest='inFastqFile',
                        required=True, help='Fastq file to trim')
    parser.add_argument('-s', '--start-pos', dest='startPos',
                        required=False, type=int, default=0,
                        help='Position of first base')
    parser.add_argument('-l', '--read-len', dest='readLen', required=True,
                        type=int, help='Length of read after trimming')
    parser.add_argument('-o', '--out-fastq', dest='outFastqFile',
                        required=True, help='Out Fastq file name')
    args = parser.parse_args()


    infile = gzip.open(args.inFastqFile, 'rt')
    outfile = gzip.open(args.outFastqFile, 'wt')

    startPos = args.startPos
    endPos = args.startPos + args.readLen
    for read in SeqIO.parse(infile, 'fastq'):
        outfile.write(read[startPos:endPos].format('fastq'))

    infile.close()
    outfile.close()

if __name__ == '__main__':
    main()
