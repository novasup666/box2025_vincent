import matplotlib.pyplot as plt
import random 
#"""

f = open("/home/novi/Documents/M1/S1/BOX/box2025_vincent/homework_01/genome_hw1.fa","r")

raw_L = f.readlines()

L = [l for l in raw_L if l[0]!=">"]

genome = "".join(L)
#"""

#genome = "azeraazterzatetazetazretareztraetraeztazetraeztraeztrazetrerztazeeeeeeeterzteraezazteraztetraezatzetraeztrazetraze"
K = 30

subwords = {i:set() for i in range(1,K+1)}

N = len(genome)

for i in range(N):
    for j in range(1,K+1):
        if i+j < N:
            subwords[j].add(genome[i:i+j])

result = [len(subwords[k]) for k in range(1,K+1)]

plt.plot(result)

plt.show()

rand_genome="".join(random.choices("ATCG",k=N))
