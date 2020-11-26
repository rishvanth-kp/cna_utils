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
print(samples)

rule all:
  input:
    expand('insert_sz/{sample}_insert_sz.pdf', sample=samples),
    expand('cna_old/{sample}_5k.pdf', sample=samples),
    expand('cna_new/{sample}_5k.pdf', sample=samples),
    expand('cna_old/{sample}_20k.pdf', sample=samples),
    expand('cna_new/{sample}_20k.pdf', sample=samples)

rule mapReads:
  input:
    ref = config['ref'],
    r1 = config['readsDir'] + '/{sample}' + config['read1Suffix'],
    r2 = config['readsDir'] + '/{sample}' + config['read1Suffix']
  output:
    sam = 'mapped_reads/{sample}.sam',
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

# rule postProcessMaps:
#   input:
#     'mapped_reads/{sample}.sam'
#   output:
#     sortBam = 'mapped_reads/{sample}_sorted.bam',
#     rmdupBam = 'mapped_reads/{sample}_rmdup.bam',
#     rmdupFstat = 'flagstat/{sample}_rmdup_flagstat.txt',
#     uniqBam = 'mapped_reads/{sample}_unique.bam',
#     uniqFstat = 'flagstat/{sample}_unique_flagstat.txt',
#     fwdSam = 'mapped_reads/{sample}_fwd.sam',
#     fwdFstat = 'flagstat/{sample}_fwd_flagstat.txt'
#   params:
#     mapq = config['minMapq']
#   threads: config['nThreads']
#   shell:
#     'samtools sort -@ {threads} -o {output.sortBam} {input}; '
#     'samtools rmdup {output.sortBam} {output.rmdupBam}; '
#     'samtools flagstat {output.rmdupBam} > {output.rmdupFstat}; '
#     'samtools view -q {params.mapq} -@ {threads} -o {output.uniqBam} '
#     '{output.rmdupBam}; '
#     'samtools flagstat {output.uniqBam} > {output.uniqFstat}; '
#     'samtools view -f 0x40 -h -@ {threads} -o {output.fwdSam} '
#     '{output.uniqBam}; '
#     'samtools flagstat {output.fwdSam} > {output.fwdFstat}'

rule postProcessMaps:
  input:
    'mapped_reads/{sample}.sam'
  output:
    collateBam = temp('mapped_reads/{sample}_coallate.bam'),
    fixmateBam = temp('mapped_reads/{sample}_fixmate.bam'),
    sortedBam = temp('mapped_reads/{sample}_sorted.bam'),
    rmdupBam = temp('mapped_reads/{sample}_rmdup.bam'),
    rmdupFstat = 'flagstat/{sample}_rmdup_flagstat.txt',
    uniqBam = 'mapped_reads/{sample}_unique.bam',
    uniqFstat = 'flagstat/{sample}_unique_flagstat.txt',
    fwdSam = 'mapped_reads/{sample}_fwd.sam',
    fwdFstat = 'flagstat/{sample}_fwd_flagstat.txt'
  params:
    mapq = config['minMapq']
  threads: config['nThreads']
  shell:
    'samtools collate -@ {threads} -o {output.collateBam} {input}; '
    'samtools fixmate -@ {threads} -m {output.collateBam} '
    '{output.fixmateBam}; '
    'samtools sort -@ {threads} -o {output.sortedBam} {output.fixmateBam}; '
    'samtools markdup -@ {threads} -r {output.sortedBam} {output.rmdupBam}; '
    'samtools flagstat {output.rmdupBam} > {output.rmdupFstat}; '
    'samtools view -@ {threads} -q {params.mapq} -F 0x800 '
    '-o {output.uniqBam} {output.rmdupBam};'
    'samtools flagstat {output.uniqBam} > {output.uniqFstat}; '
    'samtools view -@ {threads} -f 0x40 -h -o {output.fwdSam} '
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
    counts5k = 'cna_old/{sample}_5k_bincounts.bed',
    stats5k = 'cna_old/{sample}_5k_stats.txt',
    cnaPlot5k = 'cna_old/{sample}_5k.pdf',
    counts20k = 'cna_old/{sample}_20k_bincounts.bed',
    stats20k = 'cna_old/{sample}_20k_stats.txt',
    cnaPlot20k = 'cna_old/{sample}_20k.pdf'
  params:
    sampleName5k = '{sample}_5k',
    sampleName20k = '{sample}_20k',
    chromSizes = config['chromSizes'],
    binBounds5k = config['oldBinBounds5k'],
    gc5k = config['oldGc5k'],
    binBounds20k = config['oldBinBounds20k'],
    gc20k = config['oldGc20k'],
    outDir = 'cna_old'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds5k} -o {output.counts5k} -s {output.stats5k}; '
    '{config[cbsNoBad]} {output.counts5k} {params.sampleName5k} {params.gc5k} '
    '{params.outDir}; '
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds20k} -o {output.counts20k} -s {output.stats20k}; '
    '{config[cbsNoBad]} {output.counts20k} {params.sampleName20k} {params.gc20k} '
    '{params.outDir}'

rule filterDeadzone:
  input:
    'mapped_reads/{sample}_fwd.sam'
  output:
    'mapped_reads/{sample}_gz.sam'
  log:
    'logs/{sample}_filter_dz.log'
  params:
    dzBed = config['dzBed']
  shell:
    '{config[filterDz]} -i {input} -b {params.dzBed} -o {output} -v 2> {log}'

rule newCna:
  input:
    sam = 'mapped_reads/{sample}_gz.sam',
  output:
    counts5k = 'cna_new/{sample}_5k_bincounts.bed',
    stats5k = 'cna_new/{sample}_5k_stats.txt',
    cnaPlot5k = 'cna_new/{sample}_5k.pdf',
    counts20k = 'cna_new/{sample}_20k_bincounts.bed',
    stats20k = 'cna_new/{sample}_20k_stats.txt',
    cnaPlot20k = 'cna_new/{sample}_20k.pdf'
  params:
    sampleName5k = '{sample}_5k',
    sampleName20k = '{sample}_20k',
    sampleName5kNoBad = '{sample}_5k_no_bad',
    sampleName20kNoBad = '{sample}_20k_no_bad',
    chromSizes = config['chromSizes'],
    binBounds5k = config['newBinBounds5k'],
    gc5k = config['newGc5k'],
    badBins5k = config['badBins5k'],
    binBounds20k = config['newBinBounds20k'],
    gc20k = config['newGc20k'],
    badBins20k = config['badBins20k'],
    outDir = 'cna_new'
  shell:
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds5k} -o {output.counts5k} -s {output.stats5k}; '
    '{config[cbsNoBad]} {output.counts5k} {params.sampleName5k} {params.gc5k} '
    '{params.outDir}; '
    '{config[cbs]} {output.counts5k} {params.sampleName5kNoBad} {params.gc5k} '
    '{params.badBins5k} {params.outDir}; '
    '{config[binCounts]} -i {input.sam} -c {params.chromSizes} '
    '-b {params.binBounds20k} -o {output.counts20k} -s {output.stats20k}; '
    '{config[cbsNoBad]} {output.counts20k} {params.sampleName20k} {params.gc20k} '
    '{params.outDir}; '
    '{config[cbs]} {output.counts20k} {params.sampleName20kNoBad} {params.gc20k} '
    '{params.badBins20k} {params.outDir}'
