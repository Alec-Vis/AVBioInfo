""" this is the collection of functions used to map a genome and find the ORI of the genome
"""

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from main import read_fasta

    with open('e_coli-strain_973250-genome.fasta') as fp:
        name, seq = next(read_fasta(fp))
        seq = seq.replace('\n', '')

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