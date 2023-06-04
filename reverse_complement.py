#!/usr/bin/env python

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = sequence[::-1]
    reverse_complement_seq = ''.join(complement.get(base, base) for base in reverse_seq)
    return reverse_complement_seq