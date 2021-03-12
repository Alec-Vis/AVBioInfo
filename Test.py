"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np
import pandas as pd
from Lexicograph import hamming_distance


data = "AGGATAATTGATATGATTCAGTTGATATGTCTGATATGGATATGATTGATATGATGATATGCCACTGATATGTTGATATGTTGATATGGCTGATATGTAGCATTTTGATATGGTGATATGTGATATGCCTTGATATGGGGTTGATATGCTGATATGTTTTCATGATATGATGATATGTGATATGCTACTACTGATATGTGATGATATGTGTGATATGCTGATTGATATGGAGGTGATATGTTTTACTATGATATGCCCTATTTGATATGTGATATGTGATATGTGATATGTGTGCTGATATGGATTCGCTGATATGGCTGATATGTTGATATGATGATATGGATCTGTCCTATCTTGATATGCATGATATGAAATGATATGGTGATATGTCAAGATGATATGGTGATATGTGATATGCCATGATATGATGATATGGTGATATGCTGATATGTGATATGGGTGATATGGCTGTGATATGACACCCCCATCGATATGATATGTGATATGGTCCTTGATATGATGATATGTCGATGATATGTGATATGTAGATAGTGATATGTGATATGTTATGATATGCCTGATATGCAGTTGATATGTTGATATGTGATATGATTGATATGCTGATATGTGATATGCGGTCTGATATGCGCGTAGCAGTGATATGTGATATGATTATGATATGTGATATGTGATATGCTGATATGTGATATGTGGTAGCAGTGATATGGTTCGGGTGATATGTGATATGTGATATGGCTGATATGTGATATGGAATCATCTCTGTGAATGATATGACAGCTCGTGATATGTGATATGGCCCTGATATGGCTGATATGTGATATGTTGATATGCCTGATATGTATGATATGCTCCATCTGATATGTGATATGTGATATGTGATATGTGATATGTGATATGGGTCTGATATGTTGATATGTACTGATATGTGATATGGTGATATGTATGATATG"


def suffix(pat):
    """ returns the suffix of a DNA seq.
    Suffix is all characters except for the first one"""
    return pat[1:]

def base_diff(b=None):
    """ return a list of nucleotide bases that are not the current base
    if no value is given will return all possible bases (nucleotides)"""
    bases = ['a','c','g','t']
    try:
        bases.remove(b)
    except ValueError:
        pass
    return bases


def immediate_neighbors(pat):
    """ create a DNA sequence 'neighbor' 1 hamming distance away"""
    # initialize Dataframe and create first entry with column
    neighborhood = pd.DataFrame()
    neighborhood['neighbor'] = [pat]
    # iterate through the pattern seq and iteratively swap out each value with the other possible bases
    for i in range(len(pat)):
        b = pat[i]
        neighbor = list(pat)
        for base in base_diff(b):
            neighbor[i] = base
            neighborhood = neighborhood.append({'neighbor': ''.join(neighbor)}, ignore_index=True)
    return neighborhood


def neighbors(pat, d=0):
    """ recursive function that output a dataframe of all the DNA seq neighbor up to d hamming dist away"""
    if d == 0:
        return {pat}
    elif len(pat) == 1:
        return pd.DataFrame(data={'a','t','c','g'}, columns=['neighbor'])
    neighborhood = pd.DataFrame()
    suffix_neighbors = neighbors(suffix(pat), d)
    for i, string in enumerate(suffix_neighbors['neighbor']):
        #slice and concatenate dataframe
        residents = pd.DataFrame()
        if hamming_distance(suffix(pat), string) < d:
            for j in base_diff():
                # create new data frame and insert/replace into original dataframe
                new_pattern = j + string
                residents = residents.append({'neighbor': new_pattern}, ignore_index=True)
        else:
            new_pattern = pat[0] + string
            residents = residents.append({'neighbor': new_pattern}, ignore_index=True)
        # this might lead into an index out of range error
        # neighborhood = pd.concat([neighborhood.iloc[:i], residents, neighborhood.iloc[(i+1):]]).reset_index(drop=True)
        neighborhood = neighborhood.append(residents).reset_index(drop=True)
    return neighborhood



# Example:
# a = pd.DataFrame([[1, 2], [2, 3], [3, 4]])
# b = pd.DataFrame([[4, 5], [5, 6]])
# pd.concat([a.iloc[:1], b, a.iloc[(1+1):]]).reset_index(drop=True)


def iterative_neighbors(pat, d):
    pass