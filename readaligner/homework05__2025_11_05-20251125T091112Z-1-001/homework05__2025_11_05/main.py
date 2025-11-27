from readfa import readfq
from xopen import xopen
import numpy as np
import pandas as pan

key = {"A":0, "C":1, "G":2, "T":3}
comp_map = {"A":"T","C":"G","G":"C","T":"A"}


def get_sequences(file_name):
    sequences = []
    with xopen(file_name) as fasta : 
        for name,seq,_ in readfq(fasta):
            sequences.append((name,seq,len(seq)))
    return sequences

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
    matrix dynamic programming. Return distance. """
    D = np.zeros((len(x)+1, len(y)+1), dtype=int)
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

def rev_comp(seq):
    comp_seq = [comp_map[c] for c in seq]
    comp_seq.reverse()
    return "".join(comp_seq)



class mapper():

    def __init__(self,ref_seqs,k):

        self._k = k
        self._indices = []
        self._sizes = []
        self._seqs = []  
        self._names = []
        self._nOfSeqs = len(ref_seqs)
        for name,ref_seq,ref_seq_length in ref_seqs:
            self._indices.append(self._build_dic(ref_seq,ref_seq_length))
            self._sizes.append(ref_seq_length)
            self._seqs.append(ref_seq)
            self._names.append(name)


    def _build_dic (self,seq,size):
        res = {}
        k = self._k
        i = 0
        K = 4**(k-1)
        while i < size-k:
            bin_kmer = binarize(seq[i:i+k],k)
            if bin_kmer in res : 
                res[bin_kmer].append(i)
            else:
                res[bin_kmer] = [i]
            i+=k
        return res

    def seed_and_extend(self,i,read_seq,read_length):
        best = -1
        best_distance = float("inf")
        already_seen = set()
        k = self._k
        size = self._sizes[i]
        seq  = self._seqs[i]
        ref_dic = self._indices[i]
        for i in range(read_length-k+1):
            num_kmer = binarize(read_seq[i:i+k],k)
            if num_kmer in ref_dic:
                #Seeding
                indices = ref_dic[num_kmer]

                for ind in indices:
                    start_index = max(0,ind-i)
                    end_index = min(size,read_length + start_index)

                    if start_index not in already_seen:
                        score = edDistDp(read_seq,seq[start_index:end_index])
                        already_seen.add(start_index)

                        if score<best_distance:
                            best = i
                            best_distance = score
                
        return best,best_distance


    def matcher(self,read_seq,read_length,q_name):

        best = {"Query_Name": q_name,"Ref_Name":"", "Position":-1,"Strand":None, "Alignment_Score":float("inf") }
        
        for i in range(self._nOfSeqs):
            name = self._names[i]
            for rev in ["+", "-"]:
                if rev == "+":
                    index,score = self.seed_and_extend(i,read_seq,read_length)
                else:
                    index,score = self.seed_and_extend(i,rev_comp(read_seq),read_length)
                if score < best["Alignment_Score"]:
                   best = {"Query_Name": q_name,"Ref_Name": name,"Position": index, "Strand":rev, "Alignment_Score": score}
                    

            
        return best





def main(ref_fn,reads_fn,k):
    """
    Result are exported in results.tsv
    """
    res = pan.DataFrame({"Query_Name":[],"Ref_Name":[],"Position":[],"Strand":[],"Alignment_Score":[]})
    ref_seqs = get_sequences(reference_file)
    read_seqs = get_sequences(reads_file)
    m = mapper(ref_seqs,k)
    cpt = 0
    for name,seq,length in read_seqs[:100]:
        res.loc[cpt] = m.matcher(seq,length,name)
        cpt+=1
    res.to_csv('results.tsv', sep="\t")



reference_file = "ref.fa.gz"
reads_file = "reads.fq.gz"
k = 21

main(reference_file,reads_file,k)

