import fmindex
import numpy as np


def read_seq(file_name):
    f = open(file_name,"r")
    sequences = []
    l = f.readline()
    while l!="":
        if l[0]!=">":
            sequences.append(l)
        l = f.readline()

    return "".join(sequences)







a = fmindex.FMindex("BANANA")
print(a.bwt)
print(a.sa)
a.encodes_bwt_to_RLE()
print(a.bwt)
a.decodes_bwt_from_RLE()
print(a.bwt)
print(a.get_string_naive())

print("genome length: ",len(read_seq("ecoli_genome_150k.fa")))
seq = read_seq("ecoli_genome_150k.fa")
b = fmindex.FMindex(seq)
b.encodes_bwt_to_RLE()
print("RLE(BWT(genome)) length : ",len(b.bwt))
