#! /usr/bin/env python
#coding: utf-8

help_info = '''

This program can link two or more dicts with common name at column 0 together.
dict should be a tab splited file like

name_1\tattr_1\tattr_2\tattr_3
name_2\tattr_1\tattr_2\tattr_3

Usage: 

./CatAttr.py a_dict b_dict > out_dict

#the very first row will be supposed as attribute's name, passed by readline()
'''

import sys

def usage():
    sys.stderr.write(help_info)

def check_argv():
    for file_name in sys.argv[1:]:
        try:
            open(file_name).close()
        except:
            sys.stderr.write("\ncannot open %s"%file_name)
            sys.exit()

def main():
    check_argv()
    records = {}
    all_attrs = []
    content = []
    for f in sys.argv[1:]:
        handle = open(f)
        attrs = handle.readline().replace("\n","").split("\t")[1:]
        all_attrs += attrs
        for line in handle:
            if line:
                line = line.replace("\n","").split("\t")
                name = line[0]
                if name not in records:
                    records[name] = {}
                for i,attr in enumerate(attrs):
                    try:
                        assert attr not in records[name]
                    except:
                        sys.stderr.write("ERROR:name space confliction, %s appears atleast twice"%(attr))
                        sys.exit()
                    records[name][attr] = line[i+1]

    line = ["name"]+all_attrs
    content.append("\t".join(line))
    for name in records:
        line=[name]
        for attr in all_attrs:
            if attr in records[name]:
                line.append(records[name][attr])
            else:
                line.append("")
        content.append("\t".join(line))

    print "\n".join(content)

if __name__ == "__main__":
    main()
