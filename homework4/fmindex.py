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

        #public variables
        self.sa = simple_kark_sort(seq+'$')     #suffix array (str list)
        self.bwt = None                         #array of the Burrows-Wheeler transform (str np array)
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None

        self._build_bwt()

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

    def _compute_fm_count(self):
        res = Counter()
        for c in self.bwt:
            if c


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