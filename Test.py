"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np


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
