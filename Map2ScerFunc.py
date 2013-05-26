#! /usr/bin/env python
#coding: utf-8

help_info = """
This program needs name2scer_homo_file which has only two columns that map your species gene 
name to scer homo gene, the scer function comes from SGD's scer feature file.

Usage:

./Map2ScerStuff -in map_file > out_file

mapfile should looks like:

your_sp_name\tscer_name
your_sp_name\tscer_name
...

the scer_name could be system_name or function name or its alias.
"""

import sys,os
scer_feature_file = "SGD_features.tab"

def usage():
    sys.stderr.write(help_info)

def arg_verify():
    try:
        assert "-in" in sys.argv
        open(sys.argv[sys.argv.index("-in")+1]).close()
    except:
        sys.stderr.write("\nERROR: no input file\n")
        usage()
        sys.exit()

def main():
    arg_verify()
    map_file = sys.argv[sys.argv.index("-in")+1]
    sp_name2scer_name = {}
    handle = open(map_file)
    handle.readline()
    for line in handle:
        line = line.replace("\n","").split("\t")
        assert len(line) == 2
        sp_name2scer_name[line[0]] = line[1]

    all_alias2sgdid = {}
    sgdid2nameAfun = {}
    for line in open(scer_feature_file):
        line = line.replace("\n","").split("\t")
        feature_type = line[1]
        if feature_type == "ORF":
            sgdid = line[0]
            #in SGD's scer feature type file column 3,4,5 is feature name,standard name and alias
            alias = [line[3],line[4]]+line[5].split("|")
            for name in alias:
                if name:
                    all_alias2sgdid[name] = sgdid
            function = line[-1]
            sgdid2nameAfun[sgdid] = "\t".join([alias[0],alias[1],function])

    content = ["\t".join(["sp_name","scer_homo_feature_name","scer_homo_std_name","scer_homo_func"])]
    for sp_gene in sp_name2scer_name:
        if sp_name2scer_name[sp_gene] in all_alias2sgdid:
            scer_stuff = sgdid2nameAfun[all_alias2sgdid[sp_name2scer_name[sp_gene]]]
        else:
            scer_stuff = "\t\t"
        content.append("\t".join([sp_gene,scer_stuff]))
    print "\n".join(content)

if __name__ == "__main__":
    main()
