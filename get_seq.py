#! /usr/bin/env python
# coding:utf-8

'''
version 0.01 
Use YGOB database as input
version 0.02
not compatable with version 0.01
Use gff file + genome_sequence file as input
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os,sys
import argparse

parser = argparse.ArgumentParser(description='get yeast genome sequence')
parser.add_argument('genome_file',metavar='genome_fasta_file',help='reference genome file')
parser.add_argument('gff_file',metavar='gff_file',help='gff annotation file')
parser.add_argument('type_field',metavar='type_field',default='CDS',help='gene type (like CDS or gene) in gff_file')
parser.add_argument('id_field',metavar='id_field',default='ID',help='gene id (like locus_tag, Gene or ID) in gff_file')
parser.add_argument('list_file',metavar='gene_list_file',help='a gene list file')
parser.add_argument('out_put_file',metavar='out_put_file',help='out put file')
parser.add_argument('--up',dest='up_N',nargs=1,type=int,help='get up_stream N sequence of gene')
parser.add_argument('--down',dest='down_N',nargs=1,type=int,help='get down_stream N sequence of gene')
args = parser.parse_args()

class Gene():
    def __init__(self,gene_name):
        self.name = gene_name        

    def _get_seq(self,n,seq_dict,allow_overlap):
        if not allow_overlap and self.up_gene != "":
            start = max(self.start-n,1,self.up_gene.end+1)
        else:
            start = max(self.start-n,1)
        if not allow_overlap and self.down_gene != "":
            end = min(self.down_gene.start-1,self.end+n,len(seq_dict[self.seqid]))
        else:
            end = min(self.end+n,len(seq_dict[self.seqid]))
        if self.strand == "+":
            self.up_seq = seq_dict[self.seqid][start-1:self.start-1]
            self.seq = seq_dict[self.seqid][self.start-1:self.end]
            self.down_seq = seq_dict[self.seqid][self.end:end]
        if self.strand == "-":
            self.up_seq = seq_dict[self.seqid][self.end:end].reverse_complement()
            self.seq = seq_dict[self.seqid][self.start-1:self.end].reverse_complement()
            self.down_seq = seq_dict[self.seqid][start-1:self.start-1].reverse_complement()

    def up_stream(self,n,seq_dict,allow_overlap=None):
        if allow_overlap == None:
            allow_overlap = True 
        self._get_seq(n,seq_dict,allow_overlap)
        return self.up_seq

    def down_stream(self,n,seq_dict,allow_overlap=None):
        if allow_overlap == None:
            allow_overlap = True 
        self._get_seq(n,seq_dict,allow_overlap)
        return self.down_seq

    def gene_seq(self,seq_dict):
        self._get_seq(n,seq_dict,allow_overlap)
        return self.seq

def load_annotation(gff_file,type_field,id_field):
    # input gff_file
    # output genes (dict)
    # for id in genes:
    #     genes[id].name = gene_name
    #     genes[id].seqid = chromosome_num
    #     genes[id].strand = "+" or "-"
    #     genes[id].start = ORF start
    #     genes[id].end = ORF end
    #     genes[id].up_gene = up stream gene name
    #     genes[id].down_gene = down stream gene name
    #     genes[id].descrip = description

    genes = {}
    handle = open(gff_file)
    content = handle.read().split("\n")
    handle.close()
    up = ""
    for line in content:
        if line and line[0]!="#":
            seqid,source,feature,start,end,score,strand,frame,attribute = line.split("\t")
            if feature == type_field and strand == "+":
                attrs = {}
                for item in attribute.split(";"):
                    field,value = item.split("=")
                    attrs[field] = value
                name = attrs[id_field]
                genes[name] = Gene(name)
                if up != "" and genes[up].seqid == seqid:
                    genes[up].down_gene = genes[name]
                    genes[name].up_gene = genes[up]
                else:
                    genes[name].up_gene = ""
                genes[name].strand = "+"
                genes[name].start = int(start)
                genes[name].end = int(end)
                genes[name].seqid = seqid
                genes[name].descrip = attrs
                genes[name].down_gene = ""
                up = name
    up = ""
    for line in content:
        if line and line[0]!="#":
            seqid,source,feature,start,end,score,strand,frame,attribute = line.split("\t")
            if feature == type_field and strand == "-":
                attrs = {}
                for item in attribute.split(";"):
                    field,value = item.split("=")
                    attrs[field] = value
                name = attrs[id_field]
                genes[name] = Gene(name)
                if up != "" and genes[up].seqid == seqid:
                    genes[up].down_gene = genes[name]
                    genes[name].up_gene = genes[up]
                else:
                    genes[name].up_gene = ""
                genes[name].strand = "-"
                genes[name].start = int(start)
                genes[name].end = int(end)
                genes[name].seqid = seqid
                genes[name].descrip = attrs
                genes[name].down_gene = ""
                up = name
    return genes

def load_genome(genome_fasta):
    # input genome.fasta
    # return genome dict
    # genome[seqid] = seq
    genome = {}
    for record in SeqIO.parse(genome_fasta,"fasta"):
        genome[record.id] = record.seq
    return genome

def main():
    genome_records = load_genome(args.genome_file)
    gene_records = load_annotation(args.gff_file,args.type_field,args.id_field)

    gene_list = set(open(args.list_file).read().strip().split("\n"))

    records = []
    out_put = args.out_put_file
    for gene in gene_list:
        if gene in gene_records:
            if "--up" in sys.argv:
                N = args.up_N[0]
                name = "%s_%s_%d"%(gene,"up",N)
                seq = gene_records[gene].up_stream(N,genome_records)
                records.append(SeqRecord(seq,name,description=""))
            elif "--down" in sys.argv:
                N = args.down_N[0]
                name = "%s_%s_%d"%(gene,"down",N)
                seq = gene_records[gene].down_stream(N,genome_records)
                records.append(SeqRecord(seq,name,description=""))
            else:
                name = gene
                seq = gene_records[gene].gene_seq(genome_records)
                records.append(SeqRecord(seq,name,description=""))
        else:
            sys.stderr.write("%s not in Records\n"%(gene))
    SeqIO.write(records,out_put,"fasta")

if __name__ == "__main__":
    main()
