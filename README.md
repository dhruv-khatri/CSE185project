# GenomeSizer (CSE185project)

This is a project for CSE 185. GenomeSizer is a tool that takes in a .fasta file and outputs an estimated genome size. This is comparable to the kmergenie tool for genome size. 

# Install Instructions
Installation requires the ... and ... libraries to be installed
```
pip install ... 
```
Once the required libraries are installed, you can install genomesizer using the following command.
```
python setup.py install
```
If the install was successful, typing ```genomesizer --help``` should show a relevant message.

# Basic Usage
The basic usage of ```genomesizer``` is 
```
genomesizer [
```
To run ```genomesizer``` on a small test example
```
genomesizer 
```

# genomesizer options
There are 3 required inputs for the ```genomesizer``` tool. Users may specify the options below:
 - ```-k```, ```--kmer-size``` 
