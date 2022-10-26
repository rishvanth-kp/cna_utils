# Copyright (C) 2020 Kuhn-Hicks Lab, University of Southern California
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
  samples.append(sampleName)

samples = unique(samples)
print(samples)

rule all:
  input:
    expand('{}/{{sample}}_R1_001_fastqc.html'.format(config['fastqcDir']),
            sample=samples),
    expand('{}/{{sample}}_5k_seg.txt'.format(config['cnvDir']),
            sample=samples),
    expand('complexity/{sample}_lc_extrap.pdf', sample=samples)

rule fastQc:
  input:
    '{}/{}{}'.format(config['readsDir'], '{sample}',
        config['read1Suffix'])
  output:
    '{}/{{sample}}_R1_001_fastqc.html'.format(config['fastqcDir'])
  params:
    outDir = config['fastqcDir']
  shell:
    'fastqc {input} -o {params.outDir}'

rule mapReads:
  input:
    '{}/{}{}'.format(config['readsDir'], '{sample}',
      config['read1Suffix'])
  output:
    sam = temp('{}/{{sample}}.sam'.format(config['mappedDir'])),
    bam = '{}/{{sample}}.bam'.format(config['mappedDir'])
  params:
    ref = config['ref']
  threads:
    config['nThreads']
  log:
    '{}/{{sample}}_bwa.log'.format(config['logDir'])
  benchmark:
    '{}/{{sample}}_bwa.benchmark'.format(config['logDir'])
  shell:
    '{config[bwa]} mem -t {threads} {params.ref} {input} '
    '1> {output.sam} 2> {log}; '
    'samtools view -@ {threads} -b -o {output.bam} {output.sam}'

rule removePcrDups:
  input:
    '{}/{{sample}}.bam'.format(config['mappedDir'])
  output:
    sortedBam = temp('{}/{{sample}}_sorted.bam'.format(config['mappedDir'])),
    rmdupBam = '{}/{{sample}}_rmdup.bam'.format(config['mappedDir'])
  params:
    tmpDir = '{}/{{sample}}'.format(config['tmpDir'])
  threads:
    config['nThreads']
  shell:
    'samtools sort -@ {threads} -o {output.sortedBam} {input}; '
    'samtools rmdup -s {output.sortedBam} {output.rmdupBam}'

rule removeAmbig:
  input:
    '{}/{{sample}}_rmdup.bam'.format(config['mappedDir'])
  output:
    '{}/{{sample}}_unique.bam'.format(config['mappedDir'])
  params:
    minMapq = config['minMapq']
  threads:
    config['nThreads']
  shell:
    'samtools view -@ {threads} -q {params.minMapq} -F 0x800 '
    '-o {output} {input}'

rule binCounts:
  input:
    '{}/{{sample}}_unique.bam'.format(config['mappedDir'])
  output:
    counts5k = '{}/{{sample}}_5k_counts.bed'.format(config['countsDir']),
    stats5k = '{}/{{sample}}_5k_counts_stats.txt'.format(config['countsDir']),
    counts20k = '{}/{{sample}}_20k_counts.bed'.format(config['countsDir']),
    stats20k = '{}/{{sample}}_20k_counts_stats.txt'.format(config['countsDir'])
  params:
    dz = config['dzBed'],
    binBounds5k = config['binBounds5k'],
    binBounds20k = config['binBounds20k']
  shell:
    '{config[binCounts]} -i {input} -b {params.binBounds5k} '
    '-d {params.dz} -o {output.counts5k} -v > {output.stats5k}; '
    '{config[binCounts]} -i {input} -b {params.binBounds20k} '
    '-d {params.dz} -o {output.counts20k} -v > {output.stats20k}'

rule cnvProfile:
  input:
    counts5k = '{}/{{sample}}_5k_counts.bed'.format(config['countsDir']),
    counts20k = '{}/{{sample}}_20k_counts.bed'.format(config['countsDir'])
  output:
    seg5k = '{}/{{sample}}_5k_seg.txt'.format(config['cnvDir']),
    seg20k = '{}/{{sample}}_20k_seg.txt'.format(config['cnvDir'])
  params:
    gc5k = config['gc5k'],
    gc20k = config['gc20k'],
    badBins5k = config['badBins5k'],
    badBins20k = config['badBins20k'],
    outdir = '{}'.format(config['cnvDir']),
    sampleName5k = '{sample}_5k',
    sampleName20k = '{sample}_20k'
  shell:
    '{config[cnvSeg]} -b {input.counts5k} -g {params.gc5k} '
    '-e {params.badBins5k} '
    '-o {params.outdir} -n {params.sampleName5k} -v; '
    '{config[cnvSeg]} -b {input.counts20k} -g {params.gc20k} '
    '-e {params.badBins20k} '
    '-o {params.outdir} -n {params.sampleName20k} -v'

rule removeSupp:
  input:
    '{}/{{sample}}_sorted.bam'.format(config['mappedDir'])
  output:
    temp('{}/{{sample}}_nosup.bam'.format(config['mappedDir']))
  threads:
    config['nThreads']
  shell:
    'samtools view -@ {threads} -F 0x800 -o {output} {input}'

rule complexity:
  input:
    '{}/{{sample}}_nosup.bam'.format(config['mappedDir'])
  output:
    cTxt = 'complexity/{sample}_c_curve.txt',
    lcTxt = 'complexity/{sample}_lc_extrap.txt',
    lcPdf = 'complexity/{sample}_lc_extrap.pdf'
  params:
    step = config['preseqStep'],
    extrap = config['preseqExtrap']
  shell:
    '{config[preseq]} c_curve -B -P -s {params.step} -o {output.cTxt} {input}; '
    '{config[preseq]} lc_extrap -B -P -s {params.step} -e {params.extrap} '
    '-o {output.lcTxt} {input}; '
    '{config[plotLcCurve]} -c {output.cTxt} -l {output.lcTxt} '
    '-o {output.lcPdf}'  