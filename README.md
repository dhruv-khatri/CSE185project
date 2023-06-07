# GenomeSizer

Genome size is a necessary measure for analyses such as quality control, comparative genomics, and metagenomics. GenomeSizer is a tool that takes in a .fasta/.fq file and outputs the estimated genome size and a corresponding histogram file. This is comparable to the kmergenie and jellyfish tools for genome size and histogram output.

# Install Instructions
Installation requires the argparse, numpy, and matplotlib (only for histogram plotting) libraries to be installed
```
pip install argparse
pip install matplotlib
pip install numpy
```
Once the required libraries are installed, you can install genomesizer using the following command (make sure you are in desired directory).
```
git clone https://github.com/dhruv-khatri/GenomeSizer
  cd main
  make
```
If the install was successful, typing ```genomesizer.py --help``` should show a relevant message.

# Basic Usage
The basic usage of ```genomesizer``` is 
```
python genomesizer.py [inputfile] <options>
```
To run ```genomesizer``` on a small test example
```
python genomesizer.py example.fastq -k 25 -o output.histo
```

# genomesizer options
There are 3 optional inputs for the ```genomesizer``` tool. Users may specify the options below:
 - ```-k```, ```--kmer-size``` length of k-mers for size estimation
 - ```-o```, ```--output``` output file location
 - ```-c```, ```--maxcut``` maximum threshold frequency for cutting off error, determines which part of merged peaks to preserve (recommended for small datasets)
    - Value can be determined by running genomesizer once and looking at output .histo file

# File Format
There will be a .histo file output.

# Contributors
This repository was generated by Dhruv Khatri and Chris Lingunis.

Please submit a pull request with any corrections or suggestions.
