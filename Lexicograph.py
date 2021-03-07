"""
2/25/2021
@author Alec Vis

K-mer: a dna sequence of size k

This will produce a Lexicograph for a given k-mer:
    Lexicograph is a array of all the different combinations of DNA patterns for a given size k k
    the different combinations will have an index and a frequency of how often those patterns occur

The purpose of a Lexicograph is to improve the efficency of the Find_most_frequent_k-mer function

"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def load_data(path):
    ''' path is a string '''
    import pandas as pd
    data = pd.read_csv(path, delimiter='\t')
    return data


def gen_lexicograph(k):
    # itertools class for creating combinations with replacements
    from itertools import product
    # Generate a list of tuples with all the different combination of gene sequences of length k
    tuples = list(product(['a', 'c', 'g', 't'], repeat=k))  # combines the permutations with repeats
    lex_1 = np.array([''.join(tup) for tup in tuples])
    lex_1 = lex_1.reshape(len(lex_1), 1)
    freq = np.zeros(len(lex_1)).reshape(len(lex_1), 1)
    lex = np.concatenate((lex_1, freq), axis=1)
    return lex


def pattern_to_index(lex, seq):
    # find index of seq within the lexicograph
    # lexicograph will always be the first column
    index = np.where(lex == seq)
    return index[0][0]


def index_to_pattern(lex, index):
    return lex[index]


def reverse_complement(pattern):
    """takes a string and returns the complement of the string"""
    complement = ''
    for i in pattern.lower():
        if i == 'a':
            complement += 't'
        if i == 't':
            complement += 'a'
        if i == 'g':
            complement += 'c'
        if i == 'c':
            complement += 'g'
    return complement


def pattern_and_compl_match(genome, pat):
    ''' finds all patterns and complements of a pattern and return the start indexes'''
    window_size = len(pat)
    index = []
    for i in range(len(genome)):
        window = genome[i:(i + window_size)]
        if window == pat:
            index.append(i)
        elif window == reverse_complement(pat):
            index.append(i)
    return index


def fast_most_freq_kmer(lex, text):
    pass
