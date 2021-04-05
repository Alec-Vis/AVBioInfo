# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 18:26:49 2020

@author: Alec Vis

This file is a repository of the Naive pattern match functions
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import FreqPatternFinding as lx


from codetiming import Timer
t = Timer(name="class")

# Load DNA string into work space
if __name__ == '__main__':
    from main import read_fasta
    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))
        seq = seq.replace('\n', '')

def pattern_count(seq, pat):
    count = 0
    k = len(pat)
    seq = seq.lower()
    pat = pat.lower()
    for i in range(len(seq)-k):
        window = seq[i:k+i]
        if window == pat:
            count += 1
    return count


def most_frequent_kmers(seq, k):
    freq_count = set()
    count = []
    for i in range((len(seq))-k):
        window = seq[i:(i+k)]
        count.append(pattern_count(seq, window))
    max_count = max(count)
    for i in range(len(seq)-k):
        if count[i] == max_count:
            freq_count.add(seq[i:i+k])
    return freq_count

# print(f'sequence length = {len(string)}')
# print('most frequent patterns of length 13 with in the seq')
# t.start()
# print(most_frequent_kmers(string, K))
# x = round(t.stop(), 2)
# print(f'A typical bacteria genome is ~3.5 Million bases pairs therefore this would take ~{round(x*3500/60,2)} minutes for bacteria')
# print(f'a human genome is ~ 4 Billion base pairs so this would take ~{round(x*3500000/60/60,2)} hours for a human genome')

def gen_lexicograph(k):
    """ generating Lexicograph of all the different possible combinations and permuations of a string of lengh k.
    Note this is a VERY slow process and if K is larger than 10 the this process is exponentially long
    Therefore this is an inadequate for the large amount of data and string we are looking for"""
    # itertools class for creating combinations with replacements
    from itertools import product
    # Generate a list of tuples with all the different combination of gene sequences of length k
    tuples = list(product(['a', 'c', 'g', 't'], repeat=k))  # combines the permutations with repeats
    lex = [''.join(tup) for tup in tuples]
    lex_graph = pd.DataFrame()
    lex_graph['seq'] = lex
    lex_graph['freq'] = 0
    return lex_graph

def pattern_and_compl_match(genome, pat, d=0):
    """ finds all patterns and complements of a pattern and return the start indexes"""
    window_size = len(pat)
    index = []
    pat_compl = lx.reverse_complement(pat)
    for i in range(len(genome)):
        window = genome[i:(i + window_size)]
        ham_dis = lx.hamming_distance(window, pat)
        ham_dis_compl = lx.hamming_distance(window, pat_compl)
        if ham_dis <= d:
            index.append(i)
        elif ham_dis_compl <= d:
            index.append(i)
    return index

def clump_finding(genome, k, t, l):
    """input a genome and return the unique kmers that form clumps in a window of length l
    kmers form clumps when they occur at least t times within a window of length l"""
    # initialize initial window
    window = genome[0:l]
    # generate lex of initial window
    lex = lx.fast_kmer_freq_lex(window, k)
    # add clump vector to dataframe
    lex['clump'] = 0
    # search lx for freqs that are greater than t and add 1 to corresponding clump column
    lex.loc[lx.freq >= t, 'clump'] = 1
    for i in range(1, int(len(genome) - l)):
        # find first pattern and remove its freq from the lex table
        first_pat = genome[i - 1:i + k - 1]
        lex.loc[lx.seq == first_pat, 'freq'] -= 1
        # find last pattern and add one to it
        last_pat = genome[i - k:i + l]
        lex.loc[lx.seq == last_pat, 'freq'] += 1
        # increase clump if new seq had a freq >= t
        lex.loc[(lx.seq == last_pat) & (lx.freq >= t), 'clump'] += 1
    freq_pat = list(lx.loc[lx.clump >= 1, 'seq'])
    return freq_pat





















        