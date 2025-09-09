import matplotlib.pyplot as plt
import random 
#"""

f = open("/home/novi/Documents/M1/S1/BOX/box2025_vincent/homework_01/genome_hw1.fa","r")

raw_L = f.readlines()

L = [l for l in raw_L if l[0]!=">"]

genome = "".join(L)
#"""

#genome = "azeraazterzatetazetazretareztraetraeztazetraeztraeztrazetrerztazeeeeeeeterzteraezazteraztetraezatzetraeztrazetraze"
def sw_complexity(genome,N,K):
    subwords = {i:set() for i in range(1,K+1)}
    for i in range(N):
        for j in range(1,K+1):
            if i+j < N:
                subwords[j].add(genome[i:i+j])

    return [len(subwords[k]) for k in range(1,K+1)]

K = 30

N = len(genome)

rand_genome="".join(random.choices("ATCG",k=N))

def fibword(n):
    if n == 1:
        return "A"
    mn = 1
    mnp1 = 2
    w_n = "A"
    w_np1 = "AB"
    while mnp1 <n:
        tmp = w_np1
        w_np1 = w_np1 + w_n
        w_n = tmp

        tmp = mnp1
        mnp1 += mn
        mn = tmp
        
    return w_np1

q3 = sw_complexity(genome,N,K)
q4 = sw_complexity(rand_genome,N,K)
q5 = sw_complexity(fibword(N),N,K)


fig, (ax1, ax2,ax3) = plt.subplots(1, 3)
ax1.plot(q3)
ax1.set_title("subword complexity of the genome")
ax2.plot(q4)
ax2.set_title("subword complexity of a random genome")
ax3.plot(q5)
ax3.set_title("subword complexity of a fibonacci word")

plt.show()