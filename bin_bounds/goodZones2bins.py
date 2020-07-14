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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-bed', dest='inBed',
                        required=True, help='good zones in bed format')
    parser.add_argument('-b', '--bins', dest='nBins', type=int,
                        required=True, help='number of bins')
    parser.add_argument('-o', '--out-file', dest='outFile',
                        required=True, help='out file in BED format')
    args = parser.parse_args()

    # chrInfo['chrName'] = [# of mappable locations, # of bins]
    chrInfo = {'chr{}'.format(i): [0, 0] for i in range(1, 22+1)}
    chrInfo['chrX'] = [0, 0]
    chrInfo['chrY'] = [0, 0]
    print(chrInfo)

    # calculate the number of mappable locations in each chr
    totalMappable = 0
    for line in open(args.inBed):
        ch, start, end = line.strip().split()
        start = int(start)
        end = int(end)
        if ch in chrInfo:
            totalMappable += (end - start)
            chrInfo[ch][0] += (end - start)

    print(totalMappable)
    print(chrInfo)

    # allocate bins for each chr
    for k, v in chrInfo.items():
        v[1] = (v[0] / totalMappable) * args.nBins

    for k, v in chrInfo.items():
        print(k, v)

if __name__ == '__main__':
    main()
