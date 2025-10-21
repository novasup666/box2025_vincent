from pysuffixarray.core import SuffixArray
import numpy as np

def BW_transorm(seq):
    n = len(seq)
    sa = SuffixArray(seq)
    suffix_array = sa.suffix_array()
    bwt = np.char.chararray(n+1)
    for i in range(n):
        bwt[i] = seq[suffix_array[i]-1]  if suffix_array[i]>0 else '$'
    return bwt

def to_RLE(seq):
    '''
        seq in [A,T,C,G]*
    '''
    n = len(seq)
    l = []
    i = 0
    while i<n:
        current_char = seq[i]
        begining_i = i
        while i<n and seq[i] == current_char:
            i+=1
        l.append((current_char,i-begining_i))
    
    res = ".".join(map(lambda x : x[0]+str(x[1]), l))
    return res

print(BW_transorm("BANANA"))