#! /usr/bin/env python
# coding:utf-8

help_info = '''
Usage:

to get `up` `N` sequence of `species_A` genes on a `list_file`:

./get_seq.py --up N list_file > output_file

'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os,sys
import argparse

parser = argparse.ArgumentParser(description='get yeast genome sequence')
parser.add_argument('list_file',metavar='gene_list_file',help='a gene list file')
parser.add_argument('out_put_file',metavar='out_put_file',help='out put file')
parser.add_argument('--up',dest='up_N',nargs=1,type=int,help='get up_stream N sequence of gene')
parser.add_argument('--down',dest='down_N',nargs=1,type=int,help='get down_stream N sequence of gene')
args = parser.parse_args()


DATA_SET = "YGOB_dataset"
gene_list = set(open(args.list_file).read().strip().split("\n"))

class Gene():
    def __init__(self,gene_name):
        self.name = gene_name        

    def _get_seq(self,n,allow_overlap):
        if not allow_overlap and self.up_gene != "":
            start = max(self.start-n,1,self.up_gene.end+1)
        else:
            start = max(self.start-n,1)
        if not allow_overlap and self.down_gene != "":
            end = min(self.down_gene.start-1,self.end+n,len(YGOB[self.sp][self.seqid]))
        else:
            end = min(self.end+n,len(YGOB[self.sp][self.seqid]))
        if self.strand == "+":
            self.up_seq = YGOB[self.sp][self.seqid][start-1:self.start-1]
            self.seq = YGOB[self.sp][self.seqid][self.start-1:self.end]
            self.down_seq = YGOB[self.sp][self.seqid][self.end:end]
        if self.strand == "-":
            self.up_seq = YGOB[self.sp][self.seqid][self.end:end].reverse_complement()
            self.seq = YGOB[self.sp][self.seqid][self.start-1:self.end].reverse_complement()
            self.down_seq = YGOB[self.sp][self.seqid][start-1:self.start-1].reverse_complement()

    def up_stream(self,n,allow_overlap=None):
        if allow_overlap == None:
            allow_overlap = True 
        self._get_seq(n,allow_overlap)
        return self.up_seq

    def down_stream(self,n,allow_overlap=None):
        if allow_overlap == None:
            allow_overlap = True 
        self._get_seq(n,allow_overlap)
        return self.down_seq

    def gene_seq(self):
        self._get_seq(n,allow_overlap)
        return self.seq

def gene_db_init():
    # input pillars.tab and all genome.tab file
    # output gene_db
    # for record in gene_db:
    #     record.name = gene_name
    #     record.sp = it's species
    #     record.seqid = chromosome_num
    #     record.strand = "+" or "-"
    #     record.start = ORF start
    #     record.end = ORF end
    #     record.up_gene = up stream gene name
    #     record.down_gene = down stream gene name
    #     record.descrip = description

    gene_db = {}
    all_tab_files = [f for f in os.listdir(DATA_SET) if f.endswith("_genome.tab")]
    for f in all_tab_files:
        sp = f.replace("_genome.tab","")
        handle = open(DATA_SET+"/"+f)
        content = handle.read().split("\n")
        handle.close()
        up = ""
        for line in content:
            if line:
                name,strand,start,end,YGOB,seqid,short,coord,notes = line.split("\t")
                if strand == "1":
                    gene_db[name] = Gene(name)
                    gene_db[name].sp = sp
                    if up != "" and gene_db[up].seqid == int(seqid):
                        gene_db[up].down_gene = gene_db[name]
                        gene_db[name].up_gene = gene_db[up]
                    else:
                        gene_db[name].up_gene = ""
                    gene_db[name].strand = "+"
                    gene_db[name].start = int(start)
                    gene_db[name].end = int(end)
                    gene_db[name].seqid = int(seqid)
                    gene_db[name].descrip = notes.strip()
                    gene_db[name].down_gene = ""
                    up = name
        up = ""
        for line in content:
            if line:
                name,strand,start,end,YGOB,seqid,short,coord,notes = line.split("\t")
                if strand == "0":
                    gene_db[name] = Gene(name)
                    gene_db[name].sp = sp
                    if up != "" and gene_db[up].seqid == int(seqid):
                        gene_db[up].down_gene = gene_db[name]
                        gene_db[name].up_gene = gene_db[up]
                    else:
                        gene_db[name].up_gene = ""
                    gene_db[name].strand = "-"
                    gene_db[name].start = int(start)
                    gene_db[name].end = int(end)
                    gene_db[name].seqid = int(seqid)
                    gene_db[name].descrip = notes.strip()
                    gene_db[name].down_gene = ""
                    up = name
    return gene_db 

def load_seq():
    # input all sequence.fsa
    # return ygob dict
    # for sp in ygob:
    #     for seqid in ygob[sp]:
    #         ygob[sp][seqid] = seq
    ygob = {}
    genome_files = [f for f in os.listdir(DATA_SET) if f.endswith("_sequence.fsa")]
    for f in genome_files:
        sp = f.replace("_sequence.fsa","")
        ygob[sp] = {}
        fsa_file = DATA_SET+"/"+f
        for record in SeqIO.parse(fsa_file,"fasta"):
            seqid = int(record.id.split("_")[-1])
            ygob[sp][seqid] = record.seq
    return ygob

def write_homo_by_list(gene_list,n):
    records = []
    for gene in gene_list:
        record = return_record(gene,n)
        if record != "":
            records.append(record)
    file_name = "%s%s_%d_gene_%d_up.fsa"%(out_path,list_name,len(gene_list),n)
    SeqIO.write(records, file_name, "fasta")   

def main():
    records = []
    out_put = args.out_put_file
    for gene in gene_list:
        if gene in Gene_db:
            if "--up" in sys.argv:
                N = args.up_N[0]
                name = "%s_%s_%d"%(gene,"up",N)
                seq = Gene_db[gene].up_stream(N)
                records.append(SeqRecord(seq,name,description=""))
            elif "--down" in sys.argv:
                N = args.down_N[0]
                name = "%s_%s_%d"%(gene,"down",N)
                seq = Gene_db[gene].down_stream(N)
                records.append(SeqRecord(seq,name,description=""))
            else:
                name = gene
                seq = Gene_db[gene].gene_seq
                records.append(SeqRecord(seq,name,description=""))
        else:
            sys.stderr.write("%s not in Genedb\n"%(gene))
    SeqIO.write(records,out_put,"fasta")

if __name__ == "__main__":
    Gene_db = gene_db_init()
    YGOB = load_seq()
    main()
