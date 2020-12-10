![PolyA CI](https://github.com/TravisWheelerLab/polyA/workflows/PolyA%20CI/badge.svg)
![GitHub](https://img.shields.io/github/license/TravisWheelerLab/polyA)
![PyPI](https://img.shields.io/pypi/v/polyA)

# AAAAAAAAAAAAAAAA (PolyA):

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss

A tool for adjudicating between competing annotations of biological sequences.

## About

Annotation of a biological sequence is usually performed by aligning that
sequence to a database of known sequence elements. When that database contains
elements that are highly similar to each other, the proper annotation may
be ambiguous, because several entries in the database produce high-scoring
alignments. Typical annotation methods work by assigning a label based
on the candidate annotation with the highest alignment score; this can
overstate annotation certainty, mislabel boundaries, and fails to identify
large scale rearrangements or insertions within the annotated sequence.

PolyA is a software tool that adjudicates between competing alignment-based
annotations by computing estimates of annotation confidence, identifying a
trace with maximal confidence, and recursively splicing/stitching inserted
elements. PolyA communicates annotation certainty, identifies large scale
rearrangements, and detects boundaries between neighboring elements.

## Using

PolyA may be consumed as an ordinary Python package:

Install using the `pip` tool:

```
pip install polyA
```

Run from the command line:

```
polyA -h
```

or

```
python -m polyA -h
```

### Command Line

```
usage: polyA [-h] [--chunk-size CHUNK_SIZE] [--confidence]
             [--prior-counts FILE] [--shard-gap SHARD_GAP] [--sequences FILE]
             [--ultra-data FILE] [--easel-path BIN] [--ultra-path BIN]
             [--heatmap] [--log-file FILE] [--log-level LEVEL]
             [--matrix-position] [--output-path PATH] [--sequence-position]
             [--soda]
             FILE FILE

PolyA sequence adjudication tool

positional arguments:
  FILE                  Alignments file in Stockholm format
  FILE                  Substitution matrices file in PolyA matrix format

optional arguments:
  -h, --help            show this help message and exit
  --chunk-size CHUNK_SIZE
                        Size of the window in base pairs analyzed together
  --confidence          Run the confidence calculation and then exit
  --prior-counts FILE   TODO(Kaitlin)
  --shard-gap SHARD_GAP
                        Maximum alignment gap before sharding occurs
  --sequences FILE      TODO(Aubrey)
  --ultra-data FILE     TODO(Audrey)
  --easel-path BIN      Path to the esl_scorematrix program, if necessary
                        (assumed to be in PATH)
  --ultra-path BIN      Path to the ULTRA binary to use, if necessary (assumed
                        to be in PATH)
  --heatmap             Write a heatmap file to the output directory
  --log-file FILE       File to store log output in, defaults to stderr
  --log-level LEVEL     Logging level to use, 'debug' is the most noisy
  --matrix-position     Produce output in terms of the matrix position
  --output-path PATH    Directory to write output files to, defaults to
                        working directory
  --sequence-position   Produce output in terms of the target sequence
                        position
  --soda                Write a SODA visualization file to the output
                        directory
```

### Input Formats

PolyA accepts two required inputs and several optional inputs that
affect its behavior. The required inputs are an alignment file which
must contain alignments for all possible queries matching the target
sequence. This file must be in
[Stockholm](https://sonnhammer.sbc.su.se/Stockholm.html) format with
several custom metadata fields. The other required input is a set of
substitution matrices. This file uses a custom, but extremely simple
format.
	
#### Alignment File Format

Alignments for all possible queries matching the target sequence
should be contained in a single file in Stockholm format.

There are several special metadata fields that must exist for each
alignment in this file. See the example below. An explanation is
indented to the right of each field with additional detail as noted.

```
#=GF ID  MERX#DNA/TcMar-Tigger    query sequence (1)
#=GF TR  chr1:11543-28567         target sequence
#=GF SC  1153                     alignment score
#=GF SD  +                        strand
#=GF TQ  -1                       (2)
#=GF ST  127                      alignment start position on target
#=GF SP  601                      alignment stop position on target
#=GF CST 135                      alignment start position on query
#=GF CSP 628                      alignment stop position on query
#=GF FL  128                      (3)
#=GF MX  matrix_name              (4)
```

  * (1) query sequence names must be in the format 'name#family/class'
  * (2) valid values: 'q' if the alignment is on the reverse
    strand and the reversed sequence is the query; 't' if the alignment
    is on the reverse strand and the reversed sequence is the target;
    '-1' if the alignment is on the positive strand
  * (3) the flanking region of the unaligned query sequence
  * (4) the name of the substitution matrix file used to create alignment

#### Creating Alignment files from cross_match alignments

TODO(George): Ship the parsers with the PyPI package

We include a parser to convert cross_match alignment files to stockholm
alignments. 

TODO(Kaitlin): Document how the rm converter works

```
usage: python parser/cm_to_stockholm.py align_file.cm
output (can be input directly into polyA): 
    align_file.cm.sto
    align_file.cm.matrix
```

#### Substitution Matrix Files

Substitution matrix file example format (can include ambiguity codes):
* this file must include all of the matrices specified in the "#=GF MX" field of the alignment file, with corresponding and matching matrix names
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

### Output Formats

TODO(Kaitlin): This is no longer the right output format

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

esl_scorematrix is a part of the esl package in the [hmmer software suite](https://github.com/EddyRivasLab/hmmer/)

  - will compute lambda for the input score matrix
  - not needed if including lambda as a command line argument

[ULTRA](https://github.com/TravisWheelerLab/ultra)

  - will detect tandem repeat regions in a sequence FASTA file


### Using at the command line

TODO(George): Move this section up after installation
TODO(George): This help output is now incorrect

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

TODO(George): Mention flit install and the --symlink option (add a make target)

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
  - Jack Roddy
  - Travis Wheeler
