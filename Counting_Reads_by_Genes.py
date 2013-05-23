#! /usr/bin/env python
#coding: utf-8

'''
version 0.01
basic function
version 0.02
add time left
'''

import argparse
import HTSeq
import itertools
import time
import os,sys

parser = argparse.ArgumentParser(description="Use HTSeq count reads by genes",\
        usage='./%(prog)s sam_file gtf_file > out_file')
parser.add_argument('sam_file',help='input samfile')
parser.add_argument('gtf_file',help='input gtffile')

args = parser.parse_args()

sys.stderr.write("Processing..\n")
gtf_file = HTSeq.GFF_Reader(args.gtf_file)
exons = HTSeq.GenomicArrayOfSets("auto")
for feature in gtf_file:
    if feature.type == "exon":
        exons[feature.iv] += feature.name

counts = {}
for feature in gtf_file:
    if feature.type == "exon":
        counts[feature.name] = 0

cmd = "wc -l %s"%(args.sam_file)
count_line = os.popen(cmd).readline()
os.system("purge")
total_count = float(count_line.strip().split(" ")[0])
percent = 0.0

sys.stderr.write("Start Counting:\n")
sam_file = HTSeq.SAM_Reader(args.sam_file)
start_time = time.time()
for i,alnmt in enumerate(sam_file):
    if alnmt.aligned:
        intersection_set = None
        for iv2, step_set in exons[alnmt.iv].steps():
            if intersection_set is None:
                intersection_set = step_set.copy()
            else:
                intersection_set.intersection_update(step_set)
        if len(intersection_set) == 1:
            counts[list(intersection_set)[0]] += 1
    if i/total_count > percent:
        time_used = time.time() - start_time
        time_left = int(time_used/(i/total_count) - time_used)
        sys.stderr.write("%.2f %% finished\n"%(i/total_count*100))
        sys.stderr.write("%d h %d m %d s left\n"%(int(time_left/3600),\
                int((time_left%3600)/60),(time_left%3600)%60))
        percent += 0.05

for name in sorted(counts.keys()):
    print name,counts[name]
