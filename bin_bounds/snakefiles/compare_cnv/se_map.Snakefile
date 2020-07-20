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

from glob import glob
from numpy import unique

reads = glob('{}/*'.format(config['readsDir']))
samples = []
for i in reads:
  sampleName = i.replace('{}/'.format(config['readsDir']), '')
  sampleName = sampleName.replace('{}'.format(config['read1Suffix']), '')
  sampleName = sampleName.replace('{}'.format(config['read2Suffix']), '')
  samples.append(sampleName)

samples = unique(samples)

rule all:
  input:
    expand('cna_old_se/{sample}.pdf', sample=samples),
    expand('cna_new_se/{sample}.pdf', sample=samples)

rule mapReads:
  input:
    ref = config['ref'],
    r1 = config['readsDir'] + '/{sample}_R1_001.fastq.gz'
  output:
    sam = temp('mapped_reads_se/{sample}.sam'),
    flagstat = 'flagstat_se/{sample}_mapped_flagstat.txt'
  threads: config['nThreads']
  log:
    'logs/{sample}_bwa.log'
  benchmark:
    'logs/{sample}_bwa.benchmark'
  shell:
    '{config[bwa]} mem -t {threads} {input.ref} {input.r1} '
    '1> {output.sam} 2> {log}; '
    'samtools flagstat {output.sam} > {output.flagstat}'

rule postProcessMaps:
  input:
    'mapped_reads_se/{sample}.sam'
  output:
    sortBam = 'mapped_reads_se/{sample}_sorted.bam',
    rmdupBam = 'mapped_reads_se/{sample}_rmdup.bam',
    rmdupFstat = 'flagstat_se/{sample}_rmdup_flagstat.txt',
    uniqSam = 'mapped_reads_se/{sample}_unique.sam',
    uniqFstat = 'flagstat_se/{sample}_unique_flagstat.txt',
  params:
    mapq = config['minMapq']
  threads: config['nThreads']
  shell:
    'samtools sort -@ {threads} -o {output.sortBam} {input}; '
    'samtools rmdup -s {output.sortBam} {output.rmdupBam}; '
    'samtools flagstat {output.rmdupBam} > {output.rmdupFstat}; '
    'samtools view -q {params.mapq} -h -@ {threads} -o {output.uniqSam} '
    '{output.rmdupBam}; '
    'samtools flagstat {output.uniqSam} > {output.uniqFstat}; '

rule oldCna:
  input:
    sam = 'mapped_reads_se/{sample}_unique.sam',
  output:
    counts = 'cna_old_se/{sample}_bincounts.bed',
    stats = 'cna_old_se/{sample}_stats.txt',
    cnaPlot = 'cna_old_se/{sample}.pdf'
  params:
    sampleName = '{sample}',
    chromSizes = config['chromSizes'],
    binBounds = config['oldBinBounds'],
    gc = config['oldGc'],
    badBins = config['badBins'],
    outDir = 'cna_old_se'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds} -o {output.counts} -s {output.stats}; '
    '{config[cbs]} {output.counts} {params.sampleName} {params.gc} '
    '{params.badBins} {params.outDir}'

rule filterDeadzone:
  input:
    'mapped_reads_se/{sample}_unique.sam'
  output:
    'mapped_reads_se/{sample}_gz.sam'
  log:
    'logs/{sample}_filter_dz.log'
  params:
    dzBed = config['dzBed']
  shell:
    '{config[filterDz]} -i {input} -b {params.dzBed} -o {output} -v 2> {log}'

rule newCna:
  input:
    sam = 'mapped_reads_se/{sample}_gz.sam',
  output:
    counts = 'cna_new_se/{sample}_bincounts.bed',
    stats = 'cna_new_se/{sample}_stats.txt',
    cnaPlot = 'cna_new_se/{sample}.pdf'
  params:
    sampleName = '{sample}',
    chromSizes = config['chromSizes'],
    binBounds = config['newBinBounds'],
    gc = config['newGc'],
    badBins = config['badBins'],
    outDir = 'cna_new_se'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds} -o {output.counts} -s {output.stats}; '
    '{config[cbs]} {output.counts} {params.sampleName} {params.gc} '
    '{params.badBins} {params.outDir}'
