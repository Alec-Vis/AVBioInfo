# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 08:42:07 2020

@author: Alec Vis

A k-mer is a string of length k. We define Count(Text, Pattern) as the number 
of times that a k-mer Pattern appears as a substring of Text. For example,

Count(ACAACTATGCATACTATCGGGAACTATCCT,ACTAT)=3.
We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should 
account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern), our plan is to “slide a window” down Text, 
checking whether each k-mer substring of Text matches Pattern. We will therefore 
refer to the k-mer starting at position i of Text as Text(i, k). Throughout 
this book, we will often use 0-based indexing, meaning that we count starting 
at 0 instead of 1. In this case, Text begins at position 0 and ends at 
position |Text| − 1 (|Text| denotes the number of symbols in Text). 
For example, if Text = GACCATACTG, then Text(4, 3) = ATA. Note that the last 
k-mer of Text begins at position |Text| − k, e.g., the last 3-mer of GACCATACTG 
starts at position 10 − 3 = 7. This discussion results in the following pseudocode 
for computing Count(Text, Pattern).

    PatternCount(Text, Pattern)
        count ← 0
        for i ← 0 to |Text| − |Pattern|
            if Text(i, |Pattern|) = Pattern
                count ← count + 1
        return count

Implement PatternCount
Given: {DNA strings}} Text and Pattern.

Return: Count(Text, Pattern).

Sample Dataset
GCGCG
GCG
Sample Output
2

"""

# Load DNA string into work space
file = r'\rosalind_ba1a.txt'
path = r'C:\Users\Alec Vis\Rosalind\Data' + file  # r is in front to read it as a path and not a string
f = open(path, 'r')
string_raw = f.read()
f.close()

Q = string_raw.find('\n')

string_length = len(string_raw)
string = string_raw[0:Q]

pattern = string_raw[(Q + 1):(string_length - 1)]
pattern_length = len(pattern)


def pattern_count(s, p):
    count = 0
    for i in range(string_length):
        s = string[i:(i + pattern_length)]
        if s == pattern:
            count += 1
    return count


print(pattern_count(string, pattern))
