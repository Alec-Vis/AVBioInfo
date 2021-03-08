"""
@Alec Vis, 2/27/2021

test - workspace
"""
import numpy as np
import pandas as pd

# data = " AGGATAATTGATATGATTCAGTTGATATGTCTGATATGGATATGATTGATATGATGATATGCCACTGATATGTTGATATGTTGATATGGCTGATATGTAGCATTTTGATATGGTGATATGTGATATGCCTTGATATGGGGTTGATATGCTGATATGTTTTCATGATATGATGATATGTGATATGCTACTACTGATATGTGATGATATGTGTGATATGCTGATTGATATGGAGGTGATATGTTTTACTATGATATGCCCTATTTGATATGTGATATGTGATATGTGATATGTGTGCTGATATGGATTCGCTGATATGGCTGATATGTTGATATGATGATATGGATCTGTCCTATCTTGATATGCATGATATGAAATGATATGGTGATATGTCAAGATGATATGGTGATATGTGATATGCCATGATATGATGATATGGTGATATGCTGATATGTGATATGGGTGATATGGCTGTGATATGACACCCCCATCGATATGATATGTGATATGGTCCTTGATATGATGATATGTCGATGATATGTGATATGTAGATAGTGATATGTGATATGTTATGATATGCCTGATATGCAGTTGATATGTTGATATGTGATATGATTGATATGCTGATATGTGATATGCGGTCTGATATGCGCGTAGCAGTGATATGTGATATGATTATGATATGTGATATGTGATATGCTGATATGTGATATGTGGTAGCAGTGATATGGTTCGGGTGATATGTGATATGTGATATGGCTGATATGTGATATGGAATCATCTCTGTGAATGATATGACAGCTCGTGATATGTGATATGGCCCTGATATGGCTGATATGTGATATGTTGATATGCCTGATATGTATGATATGCTCCATCTGATATGTGATATGTGATATGTGATATGTGATATGTGATATGGGTCTGATATGTTGATATGTACTGATATGTGATATGGTGATATGTATGATATG"

def read_fasta(fp) -> object:
    """ Generator object of fasta file and will recursively return lines
    in the fasta file, where each line will begin with '>' followed by the sequence"""
    name, seq = None, []
    for line in fp:
        if line.startswith(">"):
            if name: yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield name, ''.join(seq)


with open('e_coli-strain_973250-genome.fasta') as fp:
    name, seq = next(read_fasta(fp))


# TODO develop the skew plot to determine the ori of replication
""" the skew plot takes adavntage of the fact that DNA polymerase only replicates in the 3' to 5' (reverse half strand)
 direction.Therefore, Okazaki fragments are created in the 5' to 3' direction and the single stranded DNA needs to 
'hangout' for a period of time making it more susceptible to mutation. the nucleotide that is extremely susceptible 
is Cytosine. Thus on the forward half strand, there is a decrease in the amount of Cytosine as it mutates into Thymine.
The repair mechanisms compound this by repairing Guanine into Adenine on the reverse half strand. 

In conclusion this implies that there will be a high C:G ratio in the reverse half strand and a low C:G ratio in the 
forward half strand. So if we walk along the DNA strand keeping track of this C:G ratio and we see a change from
low to high then we know the replication origin is at that location"""
def skew(genome):
    pass

