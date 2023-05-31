import argparse
import os
import sys

# Parse command line arguments
parser=argparse.ArgumentParser(description='Genome size estimator')

parser.add_argument('input_file', help='Path to the file')

# -k is an optional argument that presets k-mer size
parser.add_argument('-k', '--kmersize', help='length of k-mers for size estimation')
# -h is an optional argument for help message
parser.add_argument('-o', '--output', help="output file location")



args=parser.parse_args()

input_file=args.input_file
if args.kmersize is not None:
    try:
        kmersize = int(args.kmersize)
    except ValueError:
        print("Value provided to --kmersize option was not a number. Exiting now.\n")
        sys.exit(1)
else:
    kmersize = 21

    
output_file_location = args.output

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
    parser.print_help()

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
number_reads=0
genome_size=None
mean_kmer_coverage=None
kmer_size=None
# For now, assumes all reads are of equal length or filled with N to be of equal length
# TODO implement for fastq containing reads of unequal length
with open(input_file, 'r') as file:
    # Get read length
    # Process the first 4 lines of file outside loop to set read length and kmer size
    header = file.readline().strip()
    sequence = file.readline().strip()
    plus = file.readline().strip()
    quality = file.readline().strip()

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
    
    # Process first sequence to count kmers
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i : i + kmer_size]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    # Begin counting kmers for rest of file
    while True:
        # Read the next 4 lines as a group
        header = file.readline().strip()
        sequence = file.readline().strip()
        plus = file.readline().strip()
        quality = file.readline().strip()

        # Check if end of file is reached
        # Modified to account for degenerate input
        if(not header or not sequence or not plus or not quality):
            break

        # increment read count by 1
        number_reads+=1

        # Process the sequence to count kmers
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i : i + kmer_size]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    # Get Mean Kmer Coverage
    
    # Get number of unique kmers
    num_kmers_unique=0.0
    for kmer in kmer_counts.keys():
        num_kmers_unique+=1
    
    num_kmers_sequenced=0.0

    # Get total number of kmers sequenced
    for kmer in kmer_counts.keys():
        num_kmers_sequenced+=kmer_counts[kmer]

    # Calculate Mean Kmer Coverage
    mean_kmer_coverage=num_kmers_sequenced/num_kmers_unique

    # Calculate Genome Size
    genome_size=int(number_reads*(read_length-kmer_size+1)/mean_kmer_coverage)
    file.close()

# Print output genome size and other relevant statistics
# TODO: Make this message look prettier and more user-friendly
print("File processing completed")
print("Genome size: "+genome_size)
print("Mean Kmer Coverage: "+mean_kmer_coverage)
print("Number of reads: "+number_reads)
print("Read length: "+read_length)
print("Size for kmers used: "+kmer_size)
print("Generating kmer distribution histogram file")

# Add kmer counts to output .histo file

# Create kmer distribution histogram
histogram = {}

# Count kmer occurrences
for count in kmer_counts.values():
    histogram[count] = histogram.get(count, 0) + 1

# Sort the histogram dictionary by key (k-mer count)
sorted_histogram = sorted(histogram.items())

# Create the histogram file
histo_file = open(output_file_location, "w")

# Write the histogram data to the file
for count, frequency in sorted_histogram:
    histo_file.write(f"{count}\t{frequency}\n")

histo_file.close()

print("File generated succesfully")


sys.exit(0)
