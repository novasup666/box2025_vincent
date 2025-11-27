from main import BloomFilter, CountingBloomFilter 
import random as rd
from readfa import readfq
import matplotlib.pyplot as plt


def get_sequences(file_name):
    sequences = []
    with open(file_name) as fasta : 
        for name,seq,_ in readfq(fasta):
            sequences.append((name,seq,len(seq)))
    return sequences


def maintest():
    a = BloomFilter(50)
    alphabet = "ATCG"
    mers_list = ["".join([alphabet[rd.randint(0,3)] for _ in range(5)]) for _ in range(100)]
    mers = set(mers_list)
    added_mers = set(mers_list[:50])
    for mer in added_mers:
        a.insert(mer)
    nb_error = 0
    print(len(added_mers), len(mers))
    for mer in mers:
        query_result = a.mem(mer)
        nb_error += not ((mer in added_mers and query_result) or not query_result)

    print(nb_error/len(mers))

def countingtest():
    k = 21
    for i in range(1,7):
        seqs = get_sequences("data/file{i}.fa")
        data = {}
        max_seq_length = max([seq_len for _,_,seq_len in seqs])
        c_bloom_filter = CountingBloomFilter(k * max_seq_length)
        kmer_set = set()
        for _,seq,seq_len in seqs:
            kmers = [seq[i:i+k] for i in range(seq_len-k)]
            for kmer in kmers:
                c_bloom_filter.insert(kmer)
            kmer_set = union(set(kmers),kmer_set)
        
        for kmer in kmer_set:
            data[c.query(kmer)]+=1
            


maintest()