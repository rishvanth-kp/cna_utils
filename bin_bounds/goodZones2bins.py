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
import operator
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-bed', dest = 'inBed',
                        required = True, help = 'good zones in bed format')
    parser.add_argument('-c', '--chrom-sizes', dest = 'chromSizes',
                        required = True, help = 'UCSF chrom.sizes')
    parser.add_argument('-b', '--bins', dest = 'nBins', type=int,
                        required = True, help = 'number of bins')
    parser.add_argument('-o', '--out-file', dest = 'outFile',
                        required = True, help = 'out file in BED format')
    parser.add_argument('-v', '--verbose', action = 'store_true',
                        dest = 'VERBOSE', required = False,
                        help = 'print progress to stderr')
    args = parser.parse_args()

    # chrInfo['chrName'] = [# of mappable locations, # of bins, chrom len]
    chrInfo = {'chr{}'.format(i): [0] for i in range(1, 22+1)}
    chrInfo['chrX'] = [0]
    chrInfo['chrY'] = [0]

    # calculate the number of mappable locations in each chr
    if (args.VERBOSE):
        print('[READING GOODZONES]', file = sys.stderr)
    totalMappable = 0
    goodZones = {i:[] for i in chrInfo.keys()}
    for line in open(args.inBed):
        ch, start, end = line.strip().split()
        start = int(start)
        end = int(end)
        if ch in chrInfo:
            totalMappable += (end - start)
            chrInfo[ch][0] += (end - start)
            goodZones[ch].append([start, end - start])

    if (args.VERBOSE):
        print('\tMappable bases: %d' % (totalMappable), file = sys.stderr)

    # allocate bins for each chr: number of bins are rounded down,
    # and the leftover bins are added to the chrs with the largest
    # reminders
    if (args.VERBOSE):
        print('[DETERMINING BINS PER CHROMOSOME]', file = sys.stderr)
    totalBins = 0
    chrFrac = {}
    for k, v in chrInfo.items():
        bins = math.floor((v[0] / totalMappable) * args.nBins)
        chrFrac[k] = ((v[0] / totalMappable) * args.nBins) - bins
        v.append(bins)
        totalBins += bins
    while (totalBins < args.nBins):
        addToChr = max(chrFrac.items(), key = operator.itemgetter(1))[0]
        chrInfo[addToChr][1] += 1
        chrFrac[addToChr] = 0
        totalBins += 1

    # read chrom sizes
    if (args.VERBOSE):
        print('[READING CHROMOSOME SIZES]', file = sys.stderr)
    for line in open(args.chromSizes):
        ch, size = line.strip().split()
        if ch in chrInfo:
            chrInfo[ch].append(int(size))

    if (args.VERBOSE):
        print('[DETERMINING BIN BOUNDARIES]', file = sys.stderr)
        print('\t%s\t%s\t%s\t%s\t%s\t%s' %
                ('chrom', 'mappableFrac', 'mapsPerBin', 'bins',
                 'bigBins', 'smallBins'),
                file = sys.stderr)
    out = open(args.outFile, 'w')
    # determin the bin boundaries
    for k, v in chrInfo.items():
        goodZones[k].sort()
        mapsPerBin = v[0]/v[1]
        nBigBins = v[0] - (v[1] * math.floor(mapsPerBin))
        nSmallBins = v[1] - nBigBins
        assert (nSmallBins*math.floor(mapsPerBin) +
                nBigBins*math.ceil(mapsPerBin)) == v[0], 'Bin length mismatch'
        if (args.VERBOSE):
            print('\t%s\t%f\t%f\t%d\t%d\t%d' %
                    (k, v[0]/v[2], mapsPerBin, v[1], nSmallBins, nBigBins),
                    file = sys.stderr)
        binStart = 0
        binEnd = 0
        offset = 0
        for i in range(v[1] - 1):
            if (nBigBins > 0):
                mapsRemaining = math.ceil(mapsPerBin)
                nBigBins -= 1
            else:
                mapsRemaining = math.floor(mapsPerBin)
                nSmallBins -= 1

            while (mapsRemaining > goodZones[k][offset][1]):
                mapsRemaining -= goodZones[k][offset][1]
                offset += 1

            binEnd =  goodZones[k][offset][0] + mapsRemaining
            goodZones[k][offset][0] += mapsRemaining
            goodZones[k][offset][1] -= mapsRemaining
            assert (binEnd - binStart) >= math.floor(mapsPerBin), (
                'Bin too small')
            assert goodZones[k][offset][1] >= 0, 'Negaive goodZone'
            print('%s\t%d\t%d' % (k, binStart, binEnd), file = out)
            binStart = binEnd
        print('%s\t%d\t%d' % (k, binStart, v[2]), file = out)
    out.close()


if __name__ == '__main__':
    main()
