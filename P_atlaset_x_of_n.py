#! /usr/bin/env python
#coding: utf-8

help_info = '''
this program calculate probablity of at least x of n particular events happend. And n is the subset of N
x is the subset of M.

Usage:

./P_atleast_x_of_n -x int -n int -M int -N int

'''

from math import factorial
from scipy.stats import hypergeom
import argparse

parser = argparse.ArgumentParser(description="this program calculate " +\
        "probability of at least x of n particular events happend. And n "+\
        "is the subset of N, x is the subset of M.")
parser.add_argument('x',metavar='x',type=int,help="x is the particular events number of n ")
parser.add_argument('n',metavar='n',type=int,help="n is the picked put events out of N")
parser.add_argument('M',metavar='M',type=int,help="M is the particular events number of N ")
parser.add_argument('N',metavar='N',type=int,help="N is the total probably space event number ")


def verify(x,n,M,N):
    try:
        assert N>M
    except:
        sys.stderr.write("N (%d) is suppose bigger than M (%d)"%(N,M))
        sys.exit()
    try:
        assert N>n
    except:
        sys.stderr.write("N (%d) is suppose bigger than n (%d)"%(N,n))
        sys.exit()
    try:
        assert n>x
    except:
        sys.stderr.write("n (%d) is suppose bigger than x (%d)"%(n,x))
        sys.exit()
    try:
        assert M>x
    except:
        sys.stderr.write("M (%d) is suppose bigger than x (%d)"%(M,x))
        sys.exit()

def main():
    args = parser.parse_args()
    x,n,M,N = args.x,args.n,args.M,args.N
    verify(x,n,M,N)
    total_probablity = 0
    for i in range(x,min(n,M)+1,1):
        probablity = hypergeom.pmf(i,N,M,n)
        total_probablity += probablity
    print total_probablity

if __name__=="__main__":
    main()

