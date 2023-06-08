#!/usr/bin/env python

def reverse_complement(sequence):
    
    """
    Function to generate the reverse complement of a DNA sequence.

    Arguments:
    - sequence: a string representing a DNA sequence

    Returns:
    - reverse_complement_seq: a string representing the reverse complement of the input sequence
    """
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = sequence[::-1]
    reverse_complement_seq = ''.join(complement.get(base, base) for base in reverse_seq)
    return reverse_complement_seq
