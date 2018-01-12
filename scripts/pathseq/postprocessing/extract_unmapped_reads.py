#!/bin/python

# Parses SAM records from PathSeq output and writes only the unmapped reads in FASTQ format
# Use extract_unmapped_reads.sh to get unmapped reads from a BAM or sharded BAM

import sys
import argparse

parser = argparse.ArgumentParser(description="Gets unmapped reads from SAM records produced by the PathSeq pipeline (provided through stdin) and writes them to stdout in FASTQ format. Example: samtools view sample.bam | python get_unmapped_reads.py > unmapped.fq")
args = parser.parse_args()

for line in sys.stdin:
  if line[0] != '@':
    tok = line.strip().split('\t')
    isUnmapped = True
    for t in tok[11:]:
      if t.startswith('YP:Z:'):
        isUnmapped = False
        break
    if isUnmapped:
      name = tok[0]
      bases = tok[9]
      quals = tok[10]
      flag = int(tok[1])
      if flag & 0x1 and flag & 0x80:
        pair_suffix = '/2'
      else:
        pair_suffix = '/1'
      sys.stdout.write('@' + name + pair_suffix + '\n' + bases + '\n+\n' + quals + '\n')
