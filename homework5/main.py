import mmh3
import numpy as np
import math

class BloomFilter():
    def __init__(self,n,epsilon=0.01,m=0,k=0):
        '''
            m is the size of the bit array
            k is the number of hash function used
        '''
        if m and k:
            self._m = m
            self._k = k
        else:
            log_of_2 = math.log(2)
            self._m = math.ceil( - (n *math.log(epsilon)) / (log_of_2*log_of_2))
            self._k = math.ceil(self._m/n*log_of_2)
            print(f"According to given parameters n:{n} and epsilon:{epsilon}. The chosen values for k and m are k={self._k}, m={self._m}")

        self._array = np.zeros(self._m, dtype = np.uint8)
        self._hash_functions = [(lambda x : mmh3.hash(x,seed = i, signed = False)%self._m) for i in range(self._k)]

    def insert(self,e):
        for i in range(self._k):
            self._array[self._hash_functions[i](e)] = 1
    
    def mem(self,e):
        res = True
        for i in range(self._k):
            res = res and (self._array[self._hash_functions[i](e)] == 1)
        return res

class CountingBloomFilter():
    def __init__(self,n,epsilon=0.01,m=0,k=0):
        '''
            m is the size of the bit array
            k is the number of hash function used
        '''
        if m and k:
            self._m = m
            self._k = k
        else:
            log_of_2 = math.log(2)
            self._m = math.ceil( - (n *math.log(epsilon)) / (log_of_2*log_of_2))
            self._k = math.ceil(self._m/n*log_of_2)
            print(f"According to given parameters n:{n} and epsilon:{epsilon}. The chosen values for k and m are k={self._k}, m={self._m}")

        self._array = np.zeros(self._m, dtype = np.uint8)
        self._hash_functions = [(lambda x : mmh3.hash(x,seed = i, signed = False)%self._m) for i in range(self._k)]

    def insert(self,e):
        for i in range(self._k):
            self._array[self._hash_functions[i](e)] += 1
    
    def Query(self,e):
        res = float("inf")
        for i in range(self._k):
            res = min (res,self._array[self._hash_functions[i](e)] ) 
        return res