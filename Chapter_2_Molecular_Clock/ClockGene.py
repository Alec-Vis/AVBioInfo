"""
    Goal, find regulatory factors such as the ones that control the circadian rhythm
"""

# import libraries
import numpy as np
import pandas as pd
import Chapter_1_DNA_Replication_Origin.FreqPatternFinding as lx


# load data
def load_ba2a_data(string):
    data = pd.read_table(string)
    k = int(data.columns[0][0])
    d = int(data.columns[0][2])
    data = data.rename(columns={data.columns[0]: 'DNA'})
    return data, k, d


# Functions
def motif_enumeration(DNA, k, d):
    patterns = list()
    for kmer in DNA:
        # collect all kmer in a string of DNA
        pass
        for kmer_prime in DNA:
            # for each kmer all neighbors that are up to d away
            # use lx.neighbors_2 function for this
            pass

if __name__ == '__main__':
    data, k, d = load_ba2a_data("rosalind_ba2a.txt")
