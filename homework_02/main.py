N = 7

file_names = [f"homework_02/downloads/file{i}.fa" for i in range(1,N)]



def ns_csl_seq(file_name):
    f = open(file_name,"r")
    nb_seq = 0
    sequences = []
    cumulative_seq_length = 0
    l = f.readline()
    while l!="":
        if l[0]==">":
            nb_seq+=1
            sequences.append([])
        else:
            cumulative_seq_length+=len(l)
            if l[-1] == "\n":
                sequences[-1].append(l[:-1])
            else:
                sequences[-1].append(l)
        l = f.readline()
    seq = ["".join(s) for s in sequences]
    return(nb_seq,cumulative_seq_length,seq)



def cano_kmers(sequences,k):
    cano_kmers = set()

    def canonize(k,kmer):
        num_kmer = 0
        num_comp = 0
        table = {"A" : 0,"C":1,"G":2,"T":3}
        comptable = {"A":3,"C":2,"G":1,"T":0}
        for i in range(k):
            num_kmer = 4*num_kmer + table[kmer[i]]
            num_comp += comptable[kmer[i]]*(4**i)
        return min(num_kmer,num_comp)

    
    for s in sequences : 
        for i in range(len(s)-k):
            kmer = s[i:i+k]
            if "N" not in kmer : 
                cano_kmers.add(canonize(k,kmer))
    return cano_kmers

K = 20

ckmer_sets = []

for name in file_names:
    (n,l,seq) = ns_csl_seq(name)
    print(f"{name}: {n} sequences of total length {l}")
    ckmer_sets.append(cano_kmers(seq,K))

jaccardMatrix  = [[0 for _ in range(1,7)] for _ in range(1,N)]

for i in range(N-1) :
    for j in range(i,N-1) :
        a = len(ckmer_sets[i].union(ckmer_sets[j]))
        b = len(ckmer_sets[i].intersection(ckmer_sets[j]))
        jaccardMatrix[i][j] = b/a

print("The pairwise Jaccard indexes are :")
for l in jaccardMatrix:
    print(l)