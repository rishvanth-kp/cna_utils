#!/usr/bin/env python

import sys


def main():

	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]

	chrominfo = fileToDictionary("/home/cnv/CNV_dataprocessing/hg19.chrom.sizes.txt", 0)
	bins = fileToArray("/home/cnv/CNV_dataprocessing/hg19.bin.boundaries.5k.bowtie.k50.sorted.txt", 0)
	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	print len(binCounts)
	print len(bins)

	counter = 0
	dups = 0
	totalReads = 0
	prevChrompos = ""
	for x in INFILE:
		arow = x.rstrip().split("\t")
		thisChrom = arow[2]
		thisChrompos = arow[3]
		if thisChrom.find("_") > -1:
			#print thisChrom
			continue
		if thisChrom == "chrM":
			#print thisChrom
			continue
		if thisChrom == "":
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue
		#if thisChrom == "X":
		#	thisChrom = "23"
		#if thisChrom == "Y":
		#	thisChrom = "24"

		totalReads += 1
		if thisChrompos == prevChrompos:
			dups += 1
			continue
			
		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
		
		counter += 1
		#if counter % 100000 == 0:
		#	print counter
		
		indexUp = len(bins) - 1
		indexDown = 0
		indexMid = int((indexUp - indexDown) / 2.0)

		#print thisChrom, thisChrompos, thisAbspos
		while True:
			#print indexDown, indexMid, indexUp
			if thisAbspos >= long(bins[indexMid][2]):
				indexDown = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexMid
			else:
				indexUp = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexDown

			if indexUp - indexDown < 2:
				break

		#print thisChrom, thisChrompos, thisAbspos, bins[indexDown], bins[indexDown+1]
		binCounts[indexDown] += 1
		prevChrompos = thisChrompos
		
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

	binCounts.sort()
	
	STATFILE.write("TotalReads\tDupsRemoved\tReadsKept\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(dups))
	STATFILE.write("\t")
	STATFILE.write(str(counter))
	STATFILE.write("\t")
	STATFILE.write(str(binCounts[len(bins)/2]))
	STATFILE.write("\n")

	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
#	input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		if rd.has_key(id):
			#rd[id].append(arow)
			print "duplicate knowngene id = " + id
			print "arow =   " + str(arow)
			print "rd[id] = " + str(rd[id])
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
