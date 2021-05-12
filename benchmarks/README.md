# Benchmarks

This directory contains various benchmarks and performance testing tools.

## Creating Test Inputs

Follow the procedure described below to create a set of test inputs for
benchmarking.

### Chromosome File

  1. Navigate to the UCSC Genome Browser
  2. Select `Downloads`, then `Genome Data`
  3. Select a species, then `Sequence Data by Chromosome`
  4. Select the FASTA file to download

### Extract Region

Extract the region of interest using `esl-sfetch`.

```
esl-sfetch --index <FASTA input file>
esl-sfetch -c <start>..<end> <FASTA input file> <chromosome #> > <FASTA output file>
```

Get the chrom fasta file:
Go to UCSC genome browser
downloads -> genome data
human -> sequence data by chromosome
Select chrom fasta file to download

Get region of interest from chrom:
% esl-sfetch --index chrom_fasta_file.fa
% esl-sfetch -c start..end chrom_fasta_file.fa chr# > chr#_ex#.fa (or whatever you choose to name this)
In the newly created file (chr#_ex#.fa), change the header from "chrom/" to "chrom:"
Move chr#_ex#.fa to a folder called new-regions
Copy this folder to the cluster
	For me, I use
	% scp -r /Users/audrey/Desktop/new-regions ashingleton@hpc.mtech.edu:~/

Go to the cluster and copy my sbatch file: /data1/um/ashingleton/CM_test_regions.sbatch
In this file, change "/data1/um/ashingleton/new-regions" to match your new-regions path
Run this to create the alignments:
% sbatch path_to_sbatch_script

Your new-regions folder in the cluster will now have the .cm files. You'll have to copy these back to your local directory and use cm_to_stockholm.py to get the final input files

