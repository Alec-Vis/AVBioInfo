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

if __name__ == '__main__':
    from main import read_fasta

    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))

# ===========================================
""" functions needed to create the frequency array """


def symbol_to_number(symbol):
    return 'acgt'.index(symbol)


def number_to_symbol(num):
    return 'acgt'[num]


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
        return number_to_symbol(num)
    prefix_index = index // 4
    r = index % 4
    symbol = number_to_symbol(r)
    prefix_pat = number_to_pattern(prefix_index, k - 1)
    return ''.join((prefix_pat, symbol))


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


# TODO edit function to look for approximate matches
# TODO add a function that finds the most frequent approximate match through sorting
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
    return skw, skw.index(min(skw))


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


def suffix(pat):
    """ returns the suffix of a DNA seq.
    Suffix is all characters except for the first one"""
    return pat[1:]


def base_diff(b=None):
    """ return a list of nucleotide bases that are not the current base
    if no value is given will return all possible bases (nucleotides)"""
    bases = ['a', 'c', 'g', 't']
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
        return pd.DataFrame(data={'a', 't', 'c', 'g'}, columns=['neighbor'])
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
        # this might lead into an index out of range error
        # neighborhood = pd.concat([neighborhood.iloc[:i], residents, neighborhood.iloc[(i+1):]]).reset_index(drop=True)
        neighborhood = neighborhood.append(residents).reset_index(drop=True)
    return neighborhood
