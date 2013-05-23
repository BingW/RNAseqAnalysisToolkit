#! /usr/bin/env python
#coding: utf-8

import argparse
import HTSeq
import itertools

parser = argparse.ArgumentParser(description="Use HTSeq count reads by genes",\
        usage='./%(prog)s sam_file gtf_file > out_file')
parser.add_argument('sam_file',help='input samfile')
parser.add_argument('gtf_file',help='input gtffile')

args = parser.parse_args()

gtf_file = HTSeq.GFF_Reader(args.gtf_file)
exons = HTSeq.GenomicArrayOfSets("auto")
for feature in gtf_file:
    if feature.type == "exon":
        exons[feature.iv] += feature.name

counts = {}
for feature in gtf_file:
    if feature.type == "exon":
        counts[feature.name] = 0

sam_file = HTSeq.SAM_Reader(args.sam_file)
for alnmt in sam_file:
    if alnmt.aligned:
        intersection_set = None
        for iv2, step_set in exons[alnmt.iv].steps():
            if intersection_set is None:
                intersection_set = step_set.copy()
            else:
                intersection_set.intersection_update(step_set)
        if len(intersection_set) == 1:
            counts[list(intersection_set)[0]] += 1

for name in sorted(counts.keys()):
    print name,counts[name]
