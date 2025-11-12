from main import BloomFilter 
import random as rd

def maintest():
    a = BloomFilter(50)
    alphabet = "ATCG"
    mers_list = ["".join([alphabet[rd.randint(0,3)] for _ in range(5)]) for _ in range(100)]
    mers = set(mers_list)
    added_mers = set(mers_list[:50])