#! /usr/bin/env python3

from xopen import xopen 
from readfa import readfq
from collections import Counter

file_path = "reads/ecoli_sample_perfect_reads_forward.fasta.gz"

K = 31

t = 1

def get_sequences(fp):
    sequences = []
    with xopen(fp) as fasta : 
        for _,seq,_ in readfq(fasta):
            sequences.append(seq)
    return sequences

"""

def build_dbg(k,t,fp):
    sequences = get_sequences(fp)
    dbg = {}
    counts = {}
    for s in sequences:
        for i in range(len(s)-K):
            kmer = s[i:i+K]
            overlap = kmer[1:]
            next_base = s[i+K]
            if overlap in dbg:
                dbg[overlap].add(next_base)
            else:
                dbg[overlap] = set([next_base])

            if kmer in counts :
                counts[kmer] +=1
            else:
                counts[kmer] = 0 
    dbg_keys =  dbg.keys()
    for kminusone_mer in dbg_keys:
        for b in dbg[kminusone_mer]:
            if counts[kminusone_mer+b] < t : 
                dbg[kminusone_mer].remove(b)

    return dbg,counts

"""

def build_dbg(k,t,fp):
    sequences = get_sequences(fp)
    kc = Counter()
    for s in sequences:
        for i in range(len(s)-k):
            kmer = s[i:i+k]
            kc[kmer]+=1

    fkc = {k:n for k,n in kc.items() if n>t}
    return fkc 



dbg = build_dbg(K,t,file_path)



print(dbg)