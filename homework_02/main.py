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

for name in file_names:
    (n,l) = ns_csl(name)
    print(f"{name}: {n} sequences of total length {l}")
