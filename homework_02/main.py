file_names = [f"homework_02/downloads/file{i}.fa" for i in range(1,7)]
def ns_csl(file_name):
    f = open(file_name,"r")
    nb_seq = 0
    cumulative_seq_length = 0
    l = f.readline()
    while l!="":
        if l[0]==">":
            nb_seq+=1
        else:
            cumulative_seq_length+=len(l)
        l = f.readline()
    return(nb_seq,cumulative_seq_length)


def cano_kmers(sequences,k):
    cano_kmers = set()
    def canonize(kmer):
        candidate = []
        for c in kmer:
            if c == "A":
                candidate.append("C")
            if c == "C":
                candidate.append("A")
            if c == "T":
                candidate.append("G")            
            if c == "G":
                candidate.append("T")
        return min(kmer,candidate)

    for s in sequences : 
        for i in range(len(s)-k):
            cano_kmers.add(canonize(s[i:i+k]))
    return list(cano_kmers)

for name in file_names:
    (n,l) = ns_csl(name)
    print(f"{name}: {n} sequences of total length {l}")

