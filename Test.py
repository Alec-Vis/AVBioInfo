"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np
import pandas as pd
import Lexicograph as lx
from codetiming import Timer
from Naive_most_freq_kmer import pattern_count, most_frequent_kmers


t = Timer(name='class')


def suffix(pat):
    """ returns the suffix of a DNA seq.
    Suffix is all characters except for the first one"""
    return pat[1:]


def hamming_distance(s1, s2):
    """ calculate the hamming distance between two strings.
    this distance represents the number of differences between two DNA strings"""
    ham_dist = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            ham_dist += 1
        else:
            pass
    return ham_dist


def base_diff(b=None):
    """ return a list of nucleotide bases that are not the current base
    if no value is given will return all possible bases (nucleotides)"""
    bases = ['A', 'C', 'G', 'T']
    try:
        bases.remove(b)
    except ValueError:
        pass
    return bases


def neighbors_2(pat, d=0):
    """ recursive function that output a dataframe of all the DNA seq neighbor up to d hamming dist away"""
    if d == 0:
        return {pat}
    elif len(pat) == 1:
        return ['A', 'C', 'G', 'T']
    # initialize an empty dataframe
    neighborhood = list()
    # call the neighbors function and assign value
    suffix_neighbors = neighbors_2(pat[1:], d)
    # print(f'suffix_neighbors:\t{suffix_neighbors}')
    for i, string in enumerate(suffix_neighbors):
        # initialize second dataframe which will contain the set of sequences to be added to neighborhood
        residents = list()
        if hamming_distance(pat[1:], string) < d:
            for ii, j in enumerate(['A', 'C', 'G', 'T']):
                # create new array and insert/replace into original dataframe
                # print(f'j:\t\t\t{j}')
                # print(f'string:\t\t{string}')
                new_pattern = j + string
                residents.append(new_pattern)
                # print(f'new_pattern:\t{new_pattern}')
                # print(f'residents:\t\t{residents}')
        else:
            new_pattern = pat[0] + string
            residents.append(new_pattern)
            # print(f'residents:\t\t{residents}')
            # print(f'new_pattern:\t{new_pattern}')
        for resident in residents:
            neighborhood.append(resident)
        # print(f'neighborhood:\t{neighborhood}')
    return neighborhood


def neighbors(pat, d=0):
    """ recursive function that output a dataframe of all the DNA seq neighbor up to d hamming dist away"""
    if d == 0:
        return {pat}
    elif len(pat) == 1:
        return pd.DataFrame(data={'A', 'T', 'C', 'G'}, columns=['neighbor'])
    # initialize an empty dataframe
    neighborhood = pd.DataFrame()
    # call the neighbors function and assign value
    suffix_neighbors = neighbors(suffix(pat), d)
    for i, string in enumerate(suffix_neighbors['neighbor']):
        # initialize second dataframe which will contain the set of sequences to be added to neighborhood
        residents = pd.DataFrame()
        if hamming_distance(suffix(pat), string) < d:
            for j in base_diff():
                # create new data frame and insert/replace into original dataframe
                new_pattern = j + string
                residents = residents.append({'neighbor': new_pattern}, ignore_index=True)
        else:
            new_pattern = pat[0] + string
            residents = residents.append({'neighbor': new_pattern}, ignore_index=True)
        neighborhood = neighborhood.append(residents).reset_index(drop=True)
    return neighborhood


if __name__ == '__main__':
    from main import read_fasta

    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))
        seq = seq.replace('\n', '')
        pat = seq[:10]
        print('Original neighbors function runtime')
        t.start()
        neighbors(pat,3)
        t.stop()
        print('neighbors_2 function runtime')
        t.start()
        neighbors_2(pat, 3)
        t.stop()