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

./CatAttr.py a_dict b_dict -main c_dict > out_dict

#take c_dict's members as main member all member not in main will be discarded
#only one -main dict
'''

'''
version 0.01
basic function cat n dicts together

version 0.02
add -main flag which only use main members and discard all other member which 
not in main dicts
'''

import sys

def usage():
    sys.stderr.write(help_info)

def check_argv():
    try:
        assert len(sys.argv) > 1
        assert sys.argv.count("-main") <= 1
    except:
        sys.stderr.write("\nERROR: no input file\n")
        usage()
        sys.exit()

    for file_name in sys.argv[1:]:
        try:
            if file_name != "-main":
                open(file_name).close()
        except:
            sys.stderr.write("\nERROR: cannot open %s"%file_name)
            usage()
            sys.exit()

def main():
    check_argv()
    records = {}
    all_attrs = []
    content = []

    for f in sys.argv[1:]:
        if f == "-main":
            continue
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

    main_members = []
    if "-main" in sys.argv:
        f = sys.argv[sys.argv.index("-main")+1]
        handle = open(f)
        handle.readline()
        for line in handle:
            if line:
                main_members.append(line.split("\t")[0])
    else:
        main_members = records.keys()

    for name in main_members:
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
