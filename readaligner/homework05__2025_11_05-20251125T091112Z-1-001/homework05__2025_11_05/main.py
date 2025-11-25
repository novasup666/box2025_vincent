from readfa import readfq
from xopen import xopen
import numpy 

key = {"A":0, "C":1, "G":2, "T":3}
rev_map = {"A":"T","C":"G","G":"C","A":"T"}


def get_sequences(file_name):
    sequences = []
    with xopen(file_name) as fasta : 
        for name,seq,_ in readfq(fasta):
            sequences.append((name,seq,len(seq)))
    return sequences

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
    matrix dynamic programming. Return distance. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)]

def binarize(mer,k):
    if k == 1:
        return key[mer]
    else:
        res = key[mer[k-1]]+4*(binarize(mer[:k-1],k-1))
        return res

def build_dic (seq,size,k):
    res = {}
    K = 4**(k-1)
    for i in range(size-k+1):
        bin_kmer = binarize(seq[i:i+k],k)
        if bin_kmer in res : 
            res[bin_kmer].append(i)
        else:
            res[bin_kmer] = [i]
    return res

def compute_matching(ref_seq,ref_seq_length,read_seq,read_length,k):
    ref_dic = build_dic(ref_seq,ref_seq_length,k) 
    best = -1
    best_distance = float("inf")
    already_seen = set()
    for i in range(read_length-k+1):
        num_kmer = binarize(read_seq[i:i+k],k)

        if num_kmer in ref_dic:
            indices = ref_dic[num_kmer]

            for ind in indices:
                start_index = max(0,ind-i)
                end_index = min(ref_seq_length,read_length + start_index)

                if start_index not in already_seen:
                    score = edDistDp(read_seq,ref_seq[start_index:end_index])
                    already_seen.add(start_index)

                    if score<best_distance:
                        best = i
                        best_distance = score
               
    return best,best_distance


def mapper(ref_seqs,read_seq, k):
    
    read_length = len(read_seq)

    best = -1
    best_score = float("inf")
    for name,ref_seq,ref_seq_length in ref_seqs:
        matching,matching_score = compute_matching(ref_seq,ref_seq_length,read_seq,read_length,k)
        if matching_score < best_score:
            best = matching
            best_score = matching_score

def main(ref_file_name,read_file_name,k):
    ref_seqs = get_sequences(ref_file_name)
    read_seq = get_sequences(read_file_name)
    return mapper(ref_seqs,read_seq,k)

reference_file = "ref.fa.gz"
reads_file = "reads.fq.gz"

ref_seqs = get_sequences(reference_file)
read_seqs = get_sequences(reads_file)
print(read_seqs[0][1])
k = 21
for i in range(10):
    mapper(ref_seqs,read_seqs[i][1],k)