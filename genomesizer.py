import argparse
import os
import sys
import numpy as np
import copy
from reverse_complement import reverse_complement
# Parse command line arguments
parser=argparse.ArgumentParser(description='Genome size estimator')

parser.add_argument('input_file', help='Path to the file')

# -k is an optional argument that presets k-mer size
parser.add_argument('-k', '--kmersize', help='length of k-mers for size estimation')
# -o is an optional argument for output directory
parser.add_argument('-o', '--output', help="output file location")

# -c is an optional argument for maximum cut threshold for error frequencies
parser.add_argument('-c', '--maxcut', help="""maximum threshold frequency for cutting off error 
frequenices (caused by sequence errors). This option is recommended for small datasets.""")

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

maxcutfreq=args.maxcut
if (not maxcutfreq==None):
    try:
        maxcutfreq=int(maxcutfreq)
        if(maxcutfreq<=0):
            print("Value provided to --maxcut option was <= 0. Range is positive integers. Using default value")
            maxcutfreq=1000
    except ValueError:
        print("Value provided to --maxcut option was not a number. Exiting now.\n")
        sys.exit(1)
else:
    maxcutfreq=1000

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

# Check if output file location exists
if(output_file_location==None):
    print("No output file location specified. Using default.\n")
    output_file_location=input_file+".histo"

# Create output histo file
# Don't check file output extensions. Support files without extensions
try:
    output_file=open(output_file_location, "w")
    output_file.close()
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
    if not (first_sequence.startswith('A') or first_sequence.startswith('T') or first_sequence.startswith('C') or first_sequence.startswith('G')):
        print("File type is invalid. .fastq and .fq types are supported\n")
        sys.exit(1)
    first_plus=file_check.readline()
    if(not first_plus.startswith('+')):
        print("File type is invalid. .fastq and .fq types are supported\n")
        sys.exit(1)

# Our kmers and counts will be stored in the dictionary kmer_counts
avg_read_length=0
read_lengths=[]
kmer_counts={}
number_reads=0
genome_size=None
mean_kmer_coverage=None
kmer_size=None
# Get list of all read lengths
with open(input_file, 'r') as read_length_file:
    while True:
        # Read the next 4 lines as a group
        header = read_length_file.readline().strip()
        sequence = read_length_file.readline().strip()
        plus = read_length_file.readline().strip()
        quality = read_length_file.readline().strip()
        
        
        # Check if end of file is reached
        # Modified to account for degenerate input
        if(not header or not sequence or not plus or not quality):
            break
        read_lengths.append(len(sequence))
    read_length_file.close()

# Begin counting kmers
with open(input_file, 'r') as file:
    # Get read length
    # Process the first 4 lines of file outside loop to set read length and kmer size
    header = file.readline().strip()
    sequence = file.readline().strip()
    plus = file.readline().strip()
    quality = file.readline().strip()



    # Placeholder default value for kmer_size
    kmer_size=np.min(read_lengths)/2
    kmer_size = int(kmer_size)
    # set k-mer length manually if it is provided by user
    if(not kmersize==None):
        if(not kmersize>np.min(read_lengths)):
            kmer_size=kmersize
        else:
            print("User provided k-mer size larger than read length. \nUsing default value.\n")
    num_complements=0
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
            # Check for reverse complement
            rev_comp=reverse_complement(kmer)
            if(rev_comp in kmer_counts.keys()):
                num_complements+=1
                kmer_counts[rev_comp] = kmer_counts.get(rev_comp,0)+1
            else:
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
       
    print(num_complements)
    # Create kmer distribution histogram
    histogram = {}

    # Count kmer occurrences
    for count in kmer_counts.values():
        histogram[count] = histogram.get(count, 0) + 1
    
    # Fill in zero frequencies
    for i in range (1, np.max(list(kmer_counts.values()))+1):
        if(i in histogram.keys()):
            continue
        else:
            histogram[i]=0

    # Sort the histogram dictionary by key (k-mer count)
    sorted_histogram = sorted(histogram.items())

    # Get Mean Kmer Coverage
    # Trim out reads that are probably sequence errors

    trimmed_sorted_histogram=dict(sorted_histogram)
    cut_freq=-1
    max_key=max(list(trimmed_sorted_histogram.keys()))
    # Make dynamic algorithm for finding cut frequency and count
    # observed kmers manually (don't know how many error reads)

    # Loop below inaccurate in small datasets
    for frequency in range(1, max_key-1):
        if(trimmed_sorted_histogram[frequency]<trimmed_sorted_histogram[frequency+1]
           and trimmed_sorted_histogram[frequency+1]<trimmed_sorted_histogram[frequency+2]):
            cut_freq=frequency
            break
    max_cut_freq=maxcutfreq
    if (max_cut_freq>max_key):
        print("Argument provided to -c option was greater than largest frequency. Reducing to fit data size.")
        while(max_cut_freq>max_key):
            max_cut_freq=int(max_cut_freq/2)
        print("New maximum cut frequency: "+max_cut_freq)
    if (cut_freq>1000 or cut_freq>max_cut_freq):
        cut_freq=max_cut_freq
    print(cut_freq)
    error_kmers=0
    if(not cut_freq==-1):
        for i in range(1, cut_freq+1):
            error_kmers+=trimmed_sorted_histogram[i]*i
            trimmed_sorted_histogram.pop(i)

    min_key=min(list(trimmed_sorted_histogram.keys()))
    max_value = max(trimmed_sorted_histogram.values())
    max_k = [key for key, value in trimmed_sorted_histogram.items() if value == max_value][0]
    mean_kmer_coverage=max_k
    
    num_kmers_obs=0
    for frequency in range(cut_freq+1,max_key+1):
        num_kmers_obs+=trimmed_sorted_histogram[frequency]*frequency

    # Calculate average read size
    avg_read_length=np.mean(read_lengths)

    # Check to see if mean kmer coverage was obtained successfully
    if(mean_kmer_coverage==0 or mean_kmer_coverage==None):
        print("""Mean kmer coverage could not be successfully obtained. 
        Genome size cannot be estimated. Setting to -1. 
        Try looking at the output histo and setting a maximum cut threshold""")
        genome_size=-1
        mean_kmer_coverage=-1
    # Calculate Genome Size
    else:
        genome_size=int(num_kmers_obs/mean_kmer_coverage)
    file.close()

# Print output genome size and other relevant statistics
print("----------------------------------------------------------------------")
print("File processing completed")
print("Genome size: "+ str(genome_size))
print("Mean Kmer Coverage: "+ str(mean_kmer_coverage))
print("Number Non-error Kmers Observed: "+ str(num_kmers_obs))
print("Number Error Kmers Observed: "+ str(error_kmers))
print("Number of reads: "+ str(number_reads))
print("Average read length: "+ str(avg_read_length))
print("Size for kmers used: "+ str(kmer_size))
print("Generating kmer distribution histogram file...")

# Add kmer counts to output .histo file

# Create the histogram file
histo_file = open(output_file_location, "w")

# Write the histogram data to the file
for count, frequency in sorted_histogram:
    histo_file.write(f"{count}\t{frequency}\n")

histo_file.close()

print("File generated succesfully")


sys.exit(0)
