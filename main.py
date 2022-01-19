"""
2/25/2021
@author Alec Vis
This Project is

This project represents the answer to the bioinoformatics question:
    Where is the replication Origin in a DNA sequence?

Background on the Replication Origin

"""
import pandas as pd



def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        if line.startswith(">"):
            if name: yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield name, ''.join(seq)


def main():
    print('Not ready Currently preparing the functions')
    if __name__ == '__main__':
        main()
        with open('Chapter_1_DNA_Replication_Origin/human_coronavirus_HK20-42.fasta') as fp:
            name, seq = next(read_fasta(fp))
            seq = seq.replace('\n','')

