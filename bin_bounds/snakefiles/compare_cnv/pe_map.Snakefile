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
    expand('insert_sz/{sample}_insert_sz.pdf', sample=samples),
    expand('cna_old/{sample}.pdf', sample=samples),
    expand('cna_new/{sample}.pdf', sample=samples)

rule mapReads:
  input:
    ref = config['ref'],
    r1 = config['readsDir'] + '/{sample}_R1_001.fastq.gz',
    r2 = config['readsDir'] + '/{sample}_R2_001.fastq.gz'
  output:
    sam = temp('mapped_reads/{sample}.sam'),
    flagstat = 'flagstat/{sample}_mapped_flagstat.txt'
  threads: config['nThreads']
  log:
    'logs/{sample}_bwa.log'
  benchmark:
    'logs/{sample}_bwa.benchmark'
  shell:
    '{config[bwa]} mem -t {threads} {input.ref} {input.r1} {input.r2} '
    '1> {output.sam} 2> {log}; '
    'samtools flagstat {output.sam} > {output.flagstat}'

rule postProcessMaps:
  input:
    'mapped_reads/{sample}.sam'
  output:
    sortBam = 'mapped_reads/{sample}_sorted.bam',
    rmdupBam = 'mapped_reads/{sample}_rmdup.bam',
    rmdupFstat = 'flagstat/{sample}_rmdup_flagstat.txt',
    uniqBam = 'mapped_reads/{sample}_unique.bam',
    uniqFstat = 'flagstat/{sample}_unique_flagstat.txt',
    fwdSam = 'mapped_reads/{sample}_fwd.sam',
    fwdFstat = 'flagstat/{sample}_fwd_flagstat.txt'
  params:
    mapq = config['minMapq']
  threads: config['nThreads']
  shell:
    'samtools sort -@ {threads} -o {output.sortBam} {input}; '
    'samtools rmdup {output.sortBam} {output.rmdupBam}; '
    'samtools flagstat {output.rmdupBam} > {output.rmdupFstat}; '
    'samtools view -q {params.mapq} -@ {threads} -o {output.uniqBam} '
    '{output.rmdupBam}; '
    'samtools flagstat {output.uniqBam} > {output.uniqFstat}; '
    'samtools view -f 0x40 -h -@ {threads} -o {output.fwdSam} '
    '{output.uniqBam}; '
    'samtools flagstat {output.fwdSam} > {output.fwdFstat}'
    
rule insertSz:
  input:
    'mapped_reads/{sample}_unique.bam'
  output:
    pdf = 'insert_sz/{sample}_insert_sz.pdf',
    txt = 'insert_sz/{sample}_insert_sz.txt'
  log:
    'logs/{sample}_picard.log'
  benchmark:
    'logs/{sample}_picard.benchmark'
  shell:
    'java -jar {config[picard]} CollectInsertSizeMetrics I={input} '
    'O={output.txt} H={output.pdf} 2> {log}'

rule oldCna:
  input:
    sam = 'mapped_reads/{sample}_fwd.sam',
  output:
    counts = 'cna_old/{sample}_bincounts.bed',
    stats = 'cna_old/{sample}_stats.txt',
    cnaPlot = 'cna_old/{sample}.pdf'
  params:
    sampleName = '{sample}',
    chromSizes = config['chromSizes'],
    binBounds = config['oldBinBounds'],
    gc = config['oldGc'],
    badBins = config['badBins'],
    outDir = 'cna_old'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds} -o {output.counts} -s {output.stats}; '
    '{config[cbs]} {output.counts} {params.sampleName} {params.gc} '
    '{params.badBins} {params.outDir}'

rule newCna:
  input:
    sam = 'mapped_reads/{sample}_fwd.sam',
  output:
    counts = 'cna_new/{sample}_bincounts.bed',
    stats = 'cna_new/{sample}_stats.txt',
    cnaPlot = 'cna_new/{sample}.pdf'
  params:
    sampleName = '{sample}',
    chromSizes = config['chromSizes'],
    binBounds = config['newBinBounds'],
    gc = config['newGc'],
    badBins = config['badBins'],
    outDir = 'cna_new'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds} -o {output.counts} -s {output.stats}; '
    '{config[cbs]} {output.counts} {params.sampleName} {params.gc} '
    '{params.badBins} {params.outDir}'
