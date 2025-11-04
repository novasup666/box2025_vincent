#! /usr/bin/env python3

from xopen import xopen 
from readfa import readfq
from collections import Counter

file_paths = [  "reads/ecoli_sample_perfect_reads_forward.fasta.gz",
                "reads/ecoli_sample_perfect_reads.fasta.gz",
                "reads/ecoli_sample_reads_01.fasta.gz",
                "reads/ecoli_sample_reads.fasta.gz"]

K = 31
example_kmer = "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"
t = 1

memtable = {}
bases = "ATCG"
complementTable = {"A" : "T","C":"G","G":"C","T":"A"}

def get_sequences(fp):
    sequences = []
    with xopen(fp) as fasta : 
        for _,seq,_ in readfq(fasta):
            sequences.append(seq)
    return sequences


def cannonicalize(kmer):
    if kmer in memtable:
        return memtable[kmer]
    opposite =[complementTable[c] for c in kmer]
    opposite.reverse()
    memtable[kmer] = min(kmer,"".join(opposite))
    return memtable[kmer]

def build_dbg(k,t,fp):
    sequences = get_sequences(fp)
    kc = Counter()
    final_kmer_set = set()
    for s in sequences:
        for i in range(len(s)-k+1):
            kmer = cannonicalize(s[i:i+k])
            kc[kmer]+=1
            if kc[kmer] > t:
                final_kmer_set.add(kmer)
    return final_kmer_set,len(kc)

def unitig_from(dbg,kmer):
    seq = kmer
    degree = 1
    char = ""
    while degree ==1:
        degree = 0
        for c in bases:
            if cannonicalize(seq[-K+1:]+c) in dbg:
                degree+=1
                char =c
        if degree ==1:
            seq = "".join([seq,char])
    degree = 1
    char = ""
    while degree ==1:
        degree = 0
        for c in bases:
            if cannonicalize(c+seq[:K-1]) in dbg:
                degree+=1
                char =c
        if degree==1:
            seq="".join([char,seq])
    return seq
'''
for fp in file_paths:
    dbg,n = build_dbg(K,t,fp)
    print("======================")
    print(fp)
    print(f"number of {K}mers:{n} with {len(dbg)} in the graph")

for fp in file_paths[2:]:
    print("======================")
    print(fp)
    for _t in range(2,4):
        dbg,n = build_dbg(K,_t,fp)
        print(f"t={_t}, number of {K}mers:{n} with {len(dbg)} in the graph")

'''

for fp in file_paths:
    print("======================")
    print(fp)
    for _t in range(4):
        dbg,n = build_dbg(K,_t,fp)
        print(f"t={_t}, size of the unitig:{len(unitig_from(dbg, example_kmer))}")
