#!/usr/bin/env python

"""kmer_index.py: A k-mer index for indexing a text."""

__author__ = "Ben Langmead"

import bisect
pattern='GGCGCGGTGGCTCACGCCTGTAAT'



def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch += 1
                if mismatch > 2:
                    break
        if mismatch <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences

print(naive_2mm(p, t))  
print("there are", len(naive_2mm(p, t)), "occurrences")




# read in the chromosome one sequence from a FASTA file
with open('Homo_Sapiens_Chrom_1.fasta', 'r') as f:
    f.readline() # skip the first line (header)
    genome = f.read().replace('\n', '') # concatenate the rest of the lines
class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
