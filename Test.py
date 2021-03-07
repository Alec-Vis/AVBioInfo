"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np
import pandas as pd

data = " AGGATAATTGATATGATTCAGTTGATATGTCTGATATGGATATGATTGATATGATGATATGCCACTGATATGTTGATATGTTGATATGGCTGATATGTAGCATTTTGATATGGTGATATGTGATATGCCTTGATATGGGGTTGATATGCTGATATGTTTTCATGATATGATGATATGTGATATGCTACTACTGATATGTGATGATATGTGTGATATGCTGATTGATATGGAGGTGATATGTTTTACTATGATATGCCCTATTTGATATGTGATATGTGATATGTGATATGTGTGCTGATATGGATTCGCTGATATGGCTGATATGTTGATATGATGATATGGATCTGTCCTATCTTGATATGCATGATATGAAATGATATGGTGATATGTCAAGATGATATGGTGATATGTGATATGCCATGATATGATGATATGGTGATATGCTGATATGTGATATGGGTGATATGGCTGTGATATGACACCCCCATCGATATGATATGTGATATGGTCCTTGATATGATGATATGTCGATGATATGTGATATGTAGATAGTGATATGTGATATGTTATGATATGCCTGATATGCAGTTGATATGTTGATATGTGATATGATTGATATGCTGATATGTGATATGCGGTCTGATATGCGCGTAGCAGTGATATGTGATATGATTATGATATGTGATATGTGATATGCTGATATGTGATATGTGGTAGCAGTGATATGGTTCGGGTGATATGTGATATGTGATATGGCTGATATGTGATATGGAATCATCTCTGTGAATGATATGACAGCTCGTGATATGTGATATGGCCCTGATATGGCTGATATGTGATATGTTGATATGCCTGATATGTATGATATGCTCCATCTGATATGTGATATGTGATATGTGATATGTGATATGTGATATGGGTCTGATATGTTGATATGTACTGATATGTGATATGGTGATATGTATGATATG"


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
    # initialize initial window and freq_pattern
    freq_pat = []
    window = genome[0:l]
    # generate lex of initial window
    lx = fast_kmer_freq_lex(window, k)
    # add clump vector to dataframe
    lx['clump'] = 0
    # search lx for freqs that are greater than t and add 1 to corresponding clump column
    lx.loc[lx.freq >= t, 'clump'] += 1
    # TODO: create moving window, update frequency of end kmer with each slide, and update clump col
    return freq_pat

