![PolyA CI](https://github.com/TravisWheelerLab/polyA/workflows/PolyA%20CI/badge.svg)

# AAAAAAAAAAAAAAAA (PolyA):
#### a tool for adjudicating between competing annotations of biological sequences 

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss.

## About

A common annotation process compares an unannotated sequence to a collection of known sequences.
Sometimes, more than one of the queries shows a significant match to the target. Current 
annotation processes choose the 'true' match based on highest alignment score. In databases
with highly similar sequences, more than one query may match with high alignment scores, this
method is no longer reliable because it is possible that either query is the true sequence, and 
it falsely implies certainty in the true match. PolyA produces confidence estimates for 
competing annotations to eliminate implying false certainty in the annotations. In addition, 
polyA can identify instances of gene conversion and homologous recombination, and find the 
exact genomic location of the switch between sequences. It can also clarify ambiguous boundaries
between neighboring elements due to homologous over extension. It can also find sequence nesting,
identifying the inserted sequence along with the original sequence that was inserted into.


## The Algorithm

Using our confidence score analysis enables a jumping HMM approach. Allowing transitions between
competing queries identifies switches between nucleotide elements. These switches are a result of 
gene conversion, homologous recombination, neighboring elements or sequence nesting.

Our algorithm has 3 parts: 

1. Confidence calculations. This gives us the probablilites each competing query is the true 
source of the target. 
2. Position specific confidence creates a less computationally expensive version of a jpHMM, 
allowing for transitions between queries. This identifies gene conversion, homologous recombination,
nested elements, and boundaries between adjacent elements.  
3. Graph Algorithm finds nested sequences. 
	
For a more detailed description view the [poster](/publications/AlgorithmPoster.pdf)
<!-- TODO: Kaitlin - add link to paper -->

## Using

### Input File Format
	
#### Alignment Files

Alignments for all possible queries matching target sequence in single stockhold format file.

Special information fields required are:

```
#=GF ID  MERX#DNA/TcMar-Tigger              * query sequence name
#=GF TR  chr1:11543-28567                   target sequence
#=GF SC  1153                               alignment score
#=GF SD  +                                  strand
#=GF TQ  -1                                 ** see below
#=GF ST  127                                alignment start position on target
#=GF SP  601                                alignment stop position on target
#=GF CST 135                                alignment start position on query
#=GF CSP 628                                alignment stop position on query
#=GF FL  128                                *** see below
#=GF MX  matrix_name                        name of substritution matrix file used to create alignment

* query sequence names must be in the format 'name#family/class'

** TQ: 'q' if alignment is on reverse strand and the reversed sequence 
is the query. 't' if alignment is on reverse strand and the reversed 
sequence is the target. '-1' if alignment is on positive strand. 

*** FL: flanking region of unaligned query sequence.
Ex1: query sequence of length 100 aligns from 1-75, FL = 25. 
Ex2: query sequence of length 100 aligns from 10-100, FL = 9. 
```

#### Creating Alignment files from cross_match alignments

We include a parser to convert cross_match alignment files to stockholm
alignments. 

```
usage: python parser/cm_to_stockholm.py align_file.cm
output (can be input directly into polyA): 
    align_file.cm.sto
    align_file.cm.matrix
```

#### Substitution Matrix Files

Substitution matrix file example format (can include ambiguity codes):
* this file must include all of the matrices specified in the "#=GF MX" field of the alignment file, with correspoding and matching matrix names
* if lambda is not included polyA will use esl_scorematrix to calculate it for all matrices

```
matrix_name lambda(optional)
  A   G   C    T    N
  8  -6  -13  -15  -1
 -2  10  -13  -13  -1
-13  -13  10  -2   -1
-15  -13  -6   8   -1
 -1  -1   -1  -1   -1
//
matrix_name2 lambda2(optional)
  A   G   C    T    N
  8  -6  -13  -15  -1
 -2  10  -13  -13  -1
-13  -13  10  -2   -1
-15  -13  -6   8   -1
 -1  -1   -1  -1   -1
//
...
```

#### Sequence File

A FASTA file of the target sequence is needed when using ULTRA.
The target sequence must be the same genomic region that was used 
to get the cross_match alignment file.

### Output file format

```
start   stop    IDnum*   query
1       362     1111	L1PREC2_3end
363     567	2345	AluJr
568     833	3579	AluYb8
834     964	1245	AluJr
965     980	6047	L1MA4A_3end
981     1497	1111	L1PREC2_3end

* Matching IDnums correspond to partial sequences that originate from 
the same ancestral sequence.
```

#### Confidence only output file format

```
query_label         confidence
LTR40a#LTR/ERVL     0.875
LTR40b#LTR/ERVL     0.052
LTR40c#LTR/ERVL     0.001
...
```

### Extensions
#### Visualizing annotations using SODA
<!-- TODO: Kaitlin - finish this once SODA is ready -->
This section is a work in progress and will be released in the coming weeks. 

#### Prior Counts Files

The command line option --priorCounts prior_counts.txt includes prior
genome counts in confidence calculations (see paper for more details)
<!-- TODO: Kaitlin - add link to paper -->

Prior counts file example format:
```
subfamily   genome_count
AluYk2      6855
LTR38	    255
L1PA7_5end  13261
...
```

#### Using ULTRA

The command line options --seqFile seq.fasta with --ultraPath ultra_path 
will include tandem repeats in the competing annotations of the sequence.
The option --ultraOutput ultra_output.txt can also be used if ULTRA was 
ran on seq.fasta prior.


### Additional software

esl_scorematrix as a part of the esl package in the hammer software suite

  - will compute lambda for the input score matrix
  - not needed if including lambda as a command line argument

ULTRA

  - will detect tandem repeat regions in a sequence FASTA file


### Using at the command line

```
usage: python -m polyA alignFile subMatrixFile
    ARGUMENTS
        --GapInit[-25]
        --getExt[-5]
        --lambda [will calculate from substitution matrix if not included]
        --segmentsize (must be odd) [31]
        --eslPath esl_path
        --confidence - output confidence for a single annoation without running the whole algorithm
        --priorCounts prior_counts.txt
        --ultraPath ultra_path
        --seqFile genomic_region.fasta
        --ultraOutput ultra_output.txt
        --viz outfile - prints output format for SODA visualization
        --heatmap outfile - prints probability file for input into heatmap
    
    OPTIONS
        --help - display help message
        --matrixpos - prints output in terms of matrix position
        --sequencepos - prints output in terms of target sequence position
```

## Development

This project uses [Pipenv](https://pipenv.pypa.io/en/latest/), which can be
installed through Homebrew for Mac users. It must be installed before the
Makefile targets or the other commands listed in this document will work.

In order to run a command which relies on the project virtual environment, such
as `python foo.py`, it is necessary to either run `pipenv shell` first, which
will put you into a shell that has the correct Python in its `PATH`, or prefix
the command with `pipenv run` (e.g. `pipenv run python foo.py`).

### Makefile

A Makefile is available, run `make` or `make help` in the project root to see
the available targets.

### Dependencies

If you prefer not to use the Makefile, or if you need to add or remove
dependencies, the following commands will allow you to manage dependencies.

```
# Get to work (from within project directory)
# This drops you into a shell inside the project virtual environment which means
# that commands that "pipenv run" may be elided from other commands.
pipenv shell

# Fetch runtime dependencies
pipenv install

# Fetch runtime and development dependencies
pipenv install --dev

# Add a runtime dependency
pipenv install <package>

# Add a development dependency
pipenv install --dev <package>
```

### Docker Image

There is a `Dockerfile` in the repo root. The image it describes is
used for running tests in CI and can be used locally for convenience.

When dependencies change the image must be rebuilt and the new
version pushed to Docker Hub. This can be done with `make container`
if you have the correct permissions. Otherwise, ask a maintainer
to do it for you.

### Unit Tests

Unit tests use [pytest](https://pytest.org/en/latest/). The linked
documentation contains examples of how to use it to both write and run tests.
The [examples](https://pytest.org/en/latest/example/index.html) are quite
extensive.

Documentation tests are also fine for simple cases.

```
# Run tests
pipenv run python -m pytest

# or
make check-fast check-slow

# or
make check
```

## License

BSD license. See `LICENSE`.

## Authors

[Wheeler Lab](http://wheelerlab.org) at the University of Montana.

  - Kaitlin Carey
  - Audrey Shingleton
  - George Lesica
  - Travis Wheeler
