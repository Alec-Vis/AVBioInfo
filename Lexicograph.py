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
from pandas import DataFrame

if __name__ == '__main__':
    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))
    data2 = " AGGATAATTGATATGATTCAGTTGATATGTCTGATATGGATATGATTGATATGATGATATGCCACTGATATGTTGATATGTTGATATGGCTGATATGTAGCATTTTGATATGGTGATATGTGATATGCCTTGATATGGGGTTGATATGCTGATATGTTTTCATGATATGATGATATGTGATATGCTACTACTGATATGTGATGATATGTGTGATATGCTGATTGATATGGAGGTGATATGTTTTACTATGATATGCCCTATTTGATATGTGATATGTGATATGTGATATGTGTGCTGATATGGATTCGCTGATATGGCTGATATGTTGATATGATGATATGGATCTGTCCTATCTTGATATGCATGATATGAAATGATATGGTGATATGTCAAGATGATATGGTGATATGTGATATGCCATGATATGATGATATGGTGATATGCTGATATGTGATATGGGTGATATGGCTGTGATATGACACCCCCATCGATATGATATGTGATATGGTCCTTGATATGATGATATGTCGATGATATGTGATATGTAGATAGTGATATGTGATATGTTATGATATGCCTGATATGCAGTTGATATGTTGATATGTGATATGATTGATATGCTGATATGTGATATGCGGTCTGATATGCGCGTAGCAGTGATATGTGATATGATTATGATATGTGATATGTGATATGCTGATATGTGATATGTGGTAGCAGTGATATGGTTCGGGTGATATGTGATATGTGATATGGCTGATATGTGATATGGAATCATCTCTGTGAATGATATGACAGCTCGTGATATGTGATATGGCCCTGATATGGCTGATATGTGATATGTTGATATGCCTGATATGTATGATATGCTCCATCTGATATGTGATATGTGATATGTGATATGTGATATGTGATATGGGTCTGATATGTTGATATGTACTGATATGTGATATGGTGATATGTATGATATG"


def gen_lexicograph(k):
    # itertools class for creating combinations with replacements
    from itertools import product
    # Generate a list of tuples with all the different combination of gene sequences of length k
    tuples = list(product(['a', 'c', 'g', 't'], repeat=k))  # combines the permutations with repeats
    lex = [''.join(tup) for tup in tuples]
    lex_graph = pd.DataFrame()
    lex_graph['seq'] = lex;
    lex_graph['freq'] = 0
    return lex_graph


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


def hamming_distance(s1, s2):
    """ calculate the hamming distance between two strings.
    this distance represents the number of differences between two DNA strings"""
    ham_dist = 0
    s1 = s1.lower()
    s2 = s2.lower()
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            ham_dist += 1
        else:
            pass
    return ham_dist


def pattern_and_compl_match(genome, pat, d=0):
    """ finds all patterns and complements of a pattern and return the start indexes"""
    window_size = len(pat)
    index = []
    pat_compl = reverse_complement(pat)
    for i in range(len(genome)):
        window = genome[i:(i + window_size)]
        ham_dis = hamming_distance(window, pat)
        ham_dis_compl = hamming_distance(window, pat_compl)
        if ham_dis <= d:
            index.append(i)
        elif ham_dis_compl <= d:
            index.append(i)
    return index


def fast_kmer_freq_lex(genome, k):
    """ searches a genome and finds the frequency of all short gene combinations of length k (kmer)
    returns a lexicograph data frame of these frequencies for each combination"""
    lx: DataFrame = gen_lexicograph(k)
    genome = genome.lower()
    for i in range(len(genome) - k):
        window = genome[i:i + k]
        """ note: Chaining assignment does not work with data frame because we would be updating 
        this is because we would be updating a copy of the dataframe and not the dataframe itself
        therefore this must be done in a single operation with the iloc function"""
        lx.loc[lx.seq == window, 'freq'] += 1
    return lx


def clump_finding(genome, k, t, l):
    """input a genome and return the unique kmers that form clumps in a window of length l
    kmers form clumps when they occur at least t times within a window of length l"""
    # initialize initial window
    window = genome[0:l]
    # generate lex of initial window
    lx = fast_kmer_freq_lex(window, k)
    # add clump vector to dataframe
    lx['clump'] = 0
    # search lx for freqs that are greater than t and add 1 to corresponding clump column
    lx.loc[lx.freq >= t, 'clump'] = 1
    for i in range(1, int(len(genome) - l)):
        # find first pattern and remove its freq from the lex table
        first_pat = genome[i - 1:i + k - 1]
        lx.loc[lx.seq == first_pat, 'freq'] -= 1
        # find last pattern and add one to it
        last_pat = genome[i - k:i + l]
        lx.loc[lx.seq == last_pat, 'freq'] += 1
        # increase clump if new seq had a freq >= t
        lx.loc[(lx.seq == last_pat) & (lx.freq >= t), 'clump'] += 1
    freq_pat = list(lx.loc[lx.clump >= 1, 'seq'])
    return freq_pat


def skew(genome):
    """ the skew plot takes adavntage of the fact that DNA polymerase only replicates in the 3' to 5' (reverse half strand)
     direction.Therefore, Okazaki fragments are created in the 5' to 3' direction and the single stranded DNA needs to
    'hangout' for a period of time making it more susceptible to mutation. the nucleotide that is extremely susceptible
    is Cytosine. Thus on the forward half strand, there is a decrease in the amount of Cytosine as it mutates into Thymine.
    The repair mechanisms compound this by repairing Guanine into Adenine on the reverse half strand.

    In conclusion this implies that there will be a high C:G ratio in the reverse half strand and a low C:G ratio in the
    forward half strand. So if we walk along the DNA strand keeping track of this C:G ratio and we see a change from
    low to high then we know the replication origin is at that location"""
    c_count = 0
    g_count = 0
    skw = []
    genome = genome.lower()
    for i in genome:
        if i == 'c':
            c_count += 1
        elif i == 'g':
            g_count += 1
        skw.append(c_count - g_count)
    return (skw, skw.index(min(skw)))


# Generate the line plot for the skew
""" note this can take a few minutes to generate
for e coli the image is saved as a .png
Note the V shape of the graph with a minimum as ~3e6 base pairs
This is our hypothesized origin of replication"""
# print('calcuating skew')
# skw = skew(seq)
# print('skew calculation completed')
# print('generating line plot object')
# sns.lineplot(x=np.arange(len(seq)), y=skw)
# print('rendering object')
# plt.show()
