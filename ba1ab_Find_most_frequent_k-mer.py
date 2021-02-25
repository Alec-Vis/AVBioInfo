# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 18:26:49 2020

@author: Alec Vis

We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text, Pattern)
 among all k-mers. For example, "ACTAT" is a most frequent 5-mer in 
 "ACAACTATGCATCACTATCGGGAACTATCCT", 
 and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4
Sample Output
CATG GCAT

"""

# Load DNA string into work space
file = r'\rosalind_ba1b.txt'
path = r'C:\Users\Alec Vis\Rosalind\Data' + file #r is in front to read it as a path and not a string
f = open(path, 'r')
string_raw = f.read()
f.close()

# Data pre-processing
Q = string_raw.find('\n')

string_length = len(string_raw)
string = string_raw[0:Q]
K = int(string_raw[(Q+1):(string_length - 1)])

def PatternCount(s,p):
    count = 0
    pat_len = len(p)
    for i in range(len(s)):
        text = s[i:(i+pat_len)]
        if text == p:
            count += 1
    return count

def FrequentWords(s,k):
    Freq_Count = []
    for i in range((len(s))-k):
        text = s[i:(i+k)]
        Freq_Count.append(PatternCount(s,text))
    
    return Freq_Count

x = FrequentWords(string,K)























        