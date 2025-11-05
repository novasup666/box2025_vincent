# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|O|X|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)
# Code by ____

from ks import simple_kark_sort
import numpy as np
from collections import Counter

# >>==========================================================================<<
# >>                              FM-index class                              <<
# >>==========================================================================<<


class FMindex:

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False):
        #private variables
        self._seq = seq         #sequence
        self._sorted_bwt = None 

        #public variables
        self.sa = simple_kark_sort(seq+'$')     #suffix array (str list)
        self.bwt = None                         #array of the Burrows-Wheeler transform (str np array)
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None

        self._build_bwt()
        self._compute_fm_count()
        self._compute_fm_rank()
        self._compute_next_smallest_letter()

        if verbose:
            print("bwt:",self.bwt)
            print("sorted_bwt", self._sorted_bwt)
            print("fm_count:",self.fm_count)
            print("fm_rank",self.fm_rank)
            print("fm_ranks",self.fm_ranks)
            print("next_smallest_letter",self.next_smallest_letter)

    # >>===[ Attribute initialisation functions ]===============================
    def _build_bwt(self):
        '''
        Computes the Borrows Wheeler Transform of seq (of size n and given suffix_array)
        Returns a numpy char array
        '''
        n = len(self._seq)
        res = np.empty(n+1,dtype=str)
        for i in range(n+1):
            res[i] = self._seq[self.sa[i]-1]  if self.sa[i]>0 else '$'

        self.bwt = "".join(res)
        res.sort()
        self._sorted_bwt = "".join(res)

    def _compute_fm_count(self):
        res = {}
        for i in range(len(self._sorted_bwt)):
            if self._sorted_bwt[i] not in res:
                res[self._sorted_bwt[i]] = i

        self.fm_count = res
        self.fm_count["~"] = len(self.bwt)-1

    def _compute_fm_rank(self):
        cpt = Counter()
        res = [0]*len(self.bwt)
        for i in range(len(self.bwt)):
            cpt[self.bwt[i]]+=1
            res[i] = cpt[self.bwt[i]]
        self.fm_rank  = res

    def _compute_next_smallest_letter(self):
        res = {}
        letters = list(set(self._sorted_bwt))
        letters.sort()
        for i in range(len(letters)) :
            if i <len(letters)-1:
                res[letters[i]] = letters[i+1]
            else:
                res[letters[i]] = "~"
        self.next_smallest_letter = res




    # >>===[ Compression functions      ]=======================================


    def encodes_bwt_to_RLE(self):
        '''
            Compute the run length encoding of self.bwt
            Assumes that the string has no numbers in it.
        '''
        seq = self.bwt
        n = len(seq)
        l = []
        i = 0
        while i<n:
            current_char = seq[i]
            begining_i = i
            while i<n and seq[i] == current_char:
                i+=1
            l.append((current_char,i-begining_i))
        
        self.bwt = "".join(map(lambda x : x[0]+chr(x[1]), l))
        

    def decodes_bwt_from_RLE(self):
        seq = self.bwt
        n = len(seq)
        res = []
        i = 0
        while i <n :
            current = seq[i]
            i+=1
            n_of_current = ord(seq[i])
            i+=1
            res.append(current*n_of_current)
        self.bwt = "".join(res)

        

    # >>===[ Pattern matching functions ]=======================================




    def get_string_naive(self):
        '''
        Computes the sequence from the BWT matrix (in a naive way)
        '''
        n = len(self._seq)
        res = np.empty(n+1,dtype=str)
        for i in range(n+1):
            res[self.sa[i]-1] = self.bwt[i]
        return "".join(res)[:-1]

    def _lf_mapping(self, i):
        return self.fm_count[self.bwt[i]]+self.fm_rank[i]-1

    def get_string(self):
        seq = ["$"]
        i = 0
        for _ in range(len(self.bwt)-1):
            #print(next_char)
            char = self.bwt[i]
            print(seq)
            seq.append(char)
            i = self._lf_mapping(i)
        return "".join(seq[::-1])