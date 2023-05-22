import argparse
import os
import sys

# Parse command line arguments
parser=argparse.ArgumentParser(description='Genome size estimator')

parser.add_argument('input_file', help='Path to the file')

# -k is an optional argument that presets k-mer size
parser.add_argument('-k', '--kmersize', help='length of k-mers for size estimation')
# -h is an optional argument for help message
parser.add_argument('-h', '--help', help='brings up this help message', action='store_true')
parser.add_argument('-o', '--output', help="output file location")



args=parser.parse_args()

input_file=args.input_file
try:
    kmersize=int(args.k)
except ValueError:
    print("Value provided to -k option was not a number. Exiting now.\n")
    sys.exit(1)
output_file_location=args.o

"""
Begin pre-processing user input. Catch any degenerate input
and exit if the file does not meet
"""

# Check if input file exists
if(input_file==None):
    print("No file was provided\n")
    sys.exit(1)

# Check if input file can be opened
try:
    with open(input_file, 'r'):
        pass
except IOError:
    print("File could not be opened. \nUsage: python genomesizer.py <options> file_name\n")
    sys.exit(1)
except FileNotFoundError:
    print("File could not be opened. \nUsage: python genomesizer.py <options> file_name\n")
    sys.exit(1)

# Check if help message flag is provided
if(args.h):
    # Print help message
    # TODO
    print()

# Check if output file location exists
if(output_file_location==None):
    print("No output file location specified. Using default.\n")
    output_file_location=input_file+".histo"

# Create output histo file
# TODO: check output file extensions
try:
    output_file=open(output_file_location, "w")
except FileNotFoundError:
    print("Output file could not be created. Invalid name or path\n")
    sys.exit(1)




# Check if input file is in .fq or .fastq format
with open(input_file, 'r') as file_check:
    
    first_header=file_check.readline()
    if(not first_header.startswith('@')):
        print("File type is invalid. .fastq and .fq types are supported\n")
        sys.exit(1)
    first_sequence=file_check.readline()
    if(not first_sequence.startswith('A') or not first_sequence.startswith('T')
       or not first_sequence.startswith('C') or not first_sequence.startswith('G')):
        print("File type is invalid. .fastq and .fq types are supported\n")
        sys.exit(1)
    first_plus=file_check.readline()
    if(not first_plus.startswith('+')):
        print("File type is invalid. .fastq and .fq types are supported\n")
        sys.exit(1)

# Our kmers and counts will be stored in the dictionary kmer_counts
read_length=0
kmer_counts={}
# For now, assumes all reads are of equal length or filled with N to be of equal length
with open(input_file, 'r') as file:
    # Get read length
    # skip first header
    file.readline()
    sequence=file.readline()
    for base in sequence:
        read_length+=1
    # read length should now be properly set

    # Placeholder default value for kmer_size
    kmer_size=read_length/2
    # set k-mer length manually if it is provided by user
    if(not kmersize==None):
        if(not kmersize>read_length):
            kmer_size=kmersize
        else:
            print("User provided k-mer size larger than read length. \nUsing default value.\n")
    # Begin counting kmers
    #TODO
    file.close()

# Add kmer counts to output .histo file
# TODO


sys.exit(0)