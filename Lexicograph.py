"""
2/25/2021
@author Alec Vis

K-mer: a dna sequence of size k

This will produce a Lexicograph for a given k-mer:
    Lexicograph is a array of all the different combinations of DNA patterns for a given size k k
    the different combinations will have an index and a frequency of how often those patterns occur

The purpose of a Lexicograph is to improve the efficency of the Find_most_frequent_k-mer function

"""


class Lexicograph:

    def GenLexicograph(self, k):
        from itertools import permutations
        return list(permutations(['a', 'c', 'g', 't'], k))

    def PatternToIndex(self, sequence):
        pass

    def IndexToPattern(self, index, k):
        pass

