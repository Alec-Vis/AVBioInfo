"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np
import pandas as pd
import Lexicograph as lx
from codetiming import Timer
from Naive_most_freq_kmer import pattern_count, most_frequent_kmers


if __name__ == '__main__':
    from main import read_fasta
    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))

t = Timer(name='class')


def iterative_neighbors(pat, d):
    pass

seq_test = seq[:100]
k = 3


print(f'sequence length = {len(seq_test)}')
print(f'most frequent patterns of length {k} with in the seq')
t.start()
print(lx.fast_kmer_freq_lex(seq_test, k))
x = round(t.stop(), 2)

print(f'sequence length = {len(seq_test)}')
print(f'most frequent patterns of length {k} with in the seq')
t.start()
print(most_frequent_kmers(seq_test, k))
x = round(t.stop(), 2)
