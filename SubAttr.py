#! /usr/bin/env python
#coding: utf-8

help_info = '''

This program can subset of attributes from a dict.
dicts should be a tab splited file like

name_1\tattr_1\tattr_2\tattr_3
name_2\tattr_1\tattr_2\tattr_3

Usage: 

./SubAttr.py in_dict range > out_dict

range could be like

[1]          #for column (1)
[0:2]        #for column (0,1) 
[0:2,3,4:7]  #for column (0,1,3,4,5,6)

'''

import sys

def usage():
    sys.stderr.write(help_info)

def check_argv():
    try:
        open(sys.argv[1]).close()
    except:
        sys.stderr.write("\nERROR: no input file or file cannot open\n")
        usage()
        sys.exit()
    try:
        "[" in sys.argv[2] and "]" in sys.argv[2]
        num_range = sys.argv[2].replace("[","").replace("]","").replace(":",",").split(",")
        num_range = map(lambda x: int(x),num_range)
        assert sorted(num_range) == num_range
    except:
        sys.stderr.write("\nERROR: range number should be an increase sequence\n")
        usage()
        sys.exit()
        
def parse_argv():
    range_str = sys.argv[2].replace("[","").replace("]","")
    ranges = []
    for sub_range in range_str.split(","):
        if ":" in sub_range:
            ranges += range(int(sub_range.split(":")[0]),int(sub_range.split(":")[1]),1)
        else:
            ranges.append(int(sub_range))
    return set(ranges)

def main():
    check_argv()
    in_dict = sys.argv[1]
    ranges = parse_argv()
    content = []
    for line in open(in_dict):
        line = line.replace("\n","").split("\t")
        wanted = [item for i,item in enumerate(line) if i in ranges]
        content.append("\t".join(wanted))
    print "\n".join(content)

if __name__ == "__main__":
    main()

