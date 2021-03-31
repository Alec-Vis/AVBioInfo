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

# ===========================================
""" functions needed to create the frequency array """


def symbol_to_number(symbol):
    num = 'ACGT'.index(symbol)
    return num


def number_to_symbol(num):
    sym = 'ACGT'[num]
    return sym


def pattern_to_number(pat):
    """ take in a string of a DNA seq and return the lexicographically ordered index of the
    How this works:
        This takes advantage of the fact that if the last symbol was taken off of the possible DNA sequences in a
        lexicograph their order would remain the same within the array. Likewise, the last symbol will always have
        the sequence of A, C, G, T. Finally once this symbol is remove the resulting (k-1)mer is repeated four times.
        Therefore, this function will recursively remove the last symbol, find the index of this symbol (0-3) and
         this to the index of the (k-1)mer which is multiplied by four because each repeats 4 times"""
    if pat == '':
        return 0
    symbol = pat[-1]
    prefix = pat[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)


def number_to_pattern(index, k):
    """ takes the index of a sequence within a Lexicograph and the length of the DNA sequence k to
    return the corresponding kmer
    This function works in a similar way as pattern_to_number but in reverse.
    Since the last symbol will always repeat A, C, G, T taking the modulus 4 of the index will yield an index for the
    first last symbol. Recursively doing this process on the resulting quotient will continue to yield the next
    symbol over."""
    if k == 1:
        return number_to_symbol(index)
    prefix_index = index // 4
    r = index % 4
    symbol = number_to_symbol(r)
    prefix_pat = number_to_pattern(prefix_index, k - 1)
    return ''.join((prefix_pat, symbol))


# ================================================
""" efficient Most frequent kmers function """


def comp_freq_array(genome, k):
    """ Takes a DNA sequence and an integer k corresponding to the length of the DNA sequence pattern.
    The output is an array 4**k long where the numbers indicate the number of time the sequence occurred,
    the position within the arrary is the index of a lexicographly ordered list of all possible DNA sequences
    of length k"""
    freq_array = np.zeros(4 ** k)
    for i in range(len(genome) - k):
        window = genome[i:i + k]
        # print(f'window: {window}\t i: {i}')
        index = pattern_to_number(window)
        freq_array[index] += 1
    return freq_array


def fast_most_freq_pats(genome, k):
    """ This function takes in a DNA sequence and a integer k, similar to comp_freq_array.
    Returns a Set of the DNA sequences that occurred most often within the DNA sequence"""
    freq_pats = set()
    freq_array = comp_freq_array(genome, k)
    max_count = max(freq_array)
    for i in range(4 ** k):
        if freq_array[i] == max_count:
            pat = number_to_pattern(i, k)
            freq_pats.add(pat)
    return freq_pats, max_count


# =============================================
""" functions needed for reverse complement checks and approximate pattern finding"""


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
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            ham_dist += 1
        else:
            pass
    return ham_dist


def suffix(pat):
    """ returns the suffix of a DNA seq.
    Suffix is all characters except for the first one"""
    return pat[1:]


def base_diff(b=None):
    """ return a list of nucleotide bases that are not the current base
    if no value is given will return all possible bases (nucleotides)"""
    bases = ['A', 'C', 'C', 'T']
    try:
        bases.remove(b)
    except ValueError:
        pass
    return bases


def immediate_neighbors(pat):
    """ create a DNA sequence 'neighbor' 1 hamming distance away
    input a string with containing only a,c,g,or t
    output a dataframe of the strings 1 hamming distance away"""
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

def iterative_neighbors(pat, d):
    """ None recursive version of the the Neighbors function
    this repeatedly calls the immediate neighbors function a number of times equal to the hamming distance away
    Although it is not as efficient as the neighbors function because the immediate neighbors will make repeats
    and these need to be dropped at each iteration of the loop"""
    neighborhood = pd.DataFrame()
    neighborhood['neighbor'] = [pat]
    for i in range(d):
        for index, row in neighborhood.iterrows():
            seq = row['neighbor']
            _neighborhood = immediate_neighbors(seq)
            neighborhood = pd.concat([neighborhood, _neighborhood])
            neighborhood = neighborhood.drop_duplicates()
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


# ===============================================
""" efficient most freq kmers function used to find approximate match and reverse complements"""


# TODO add a function that finds the most frequent approximate match through sorting
def comp_freq_array_with_mismatch(genome, k, d):
    """ searches a genome and finds the frequency of all short gene combinations of length k (kmer)
    returns a lexicograph numpy array of these frequencies for each combination"""
    freq_array = np.zeros(4 ** k)
    genome = genome.upper()
    for i in range(len(genome) - k):
        window = genome[i:i + k]
        neighborhood = neighbors(window, d)
        for index, row in neighborhood.iterrows():
            approximate_neighbor = row['neighbor']
            j = pattern_to_number(approximate_neighbor)
            freq_array[j] += 1
    return freq_array


# ==========================================
if __name__ == '__main__':
    from main import read_fasta

    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))
        seq = seq.replace('\n', '')
