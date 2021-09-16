![PolyA CI](https://github.com/TravisWheelerLab/polyA/workflows/PolyA%20CI/badge.svg)
![GitHub](https://img.shields.io/github/license/TravisWheelerLab/polyA)
![PyPI](https://img.shields.io/pypi/v/polyA)

# AAAAAAAAAAAAAAAA (PolyA):

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss

A tool for adjudicating between competing annotations of biological sequences.

## About

[Preprint on
bioRxiv](https://www.biorxiv.org/content/10.1101/2021.02.13.430877v1)

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

### Docker

A runner image can be found on
[Docker Hub](https://hub.docker.com/repository/docker/traviswheelerlab/polya).
To use it, see the README included alongside the image.

### Command Line

Command line usage is available with `polyA -h`. It is also included below for
convenience.

```
usage: polyA [-h] [-v] [--chunk-size CHUNK_SIZE] [--confidence]
             [--prior-counts FILE] [--shard-gap SHARD_GAP] [--sequences SEQS]
             [--ultra-data FILE] [--easel-path BIN] [--ultra-path BIN]
             [--log-file LOG] [--log-level LEVEL] [--matrix-position]
             [--output-path PATH] [--sequence-position] [--soda]
             ALIGNMENTS MATRICES

PolyA sequence adjudication tool

positional arguments:
  ALIGNMENTS              alignments file in Stockholm format
  MATRICES                substitution matrices file in PolyA matrix format

optional arguments:
  -h, --help              show this help message and exit
  -v, --version           show version and exit
  --chunk-size CHUNK_SIZE
                          size of the window in base pairs analyzed together
  --confidence            run the confidence calculation and then exit
  --prior-counts FILE     file containing query genomic counts
  --shard-gap SHARD_GAP
                          maximum alignment gap before sharding occurs
  --sequences SEQS        FASTA file of the target sequence for using ULTRA
  --ultra-data FILE       text file of the output from ULTRA ran on the FASTA file
                          of the target sequence
  --easel-path BIN        path to the esl_scorematrix program, if necessary
                          (assumed to be in PATH)
  --ultra-path BIN        path to the ULTRA binary to use, if necessary (assumed
                          to be in PATH)
  --log-file LOG          file to store log output in, defaults to stderr
  --log-level LEVEL       logging level to use, 'debug' is the most noisy
  --matrix-position       produce output in terms of the matrix position
  --output-path PATH      directory to write output files to, defaults to
                          working directory
  --sequence-position     produce output in terms of the target sequence
                          position
  --soda                  write a SODA visualization file to the output
                          directory
  --complexity-adjustment complexity-adjust alignment scores (scores of matches
                          between sequence regions of biassed nucleotide
                          composition are adjusted downwards)
```

### Input Formats

PolyA accepts two required inputs and several optional inputs that affect
its behavior. The required inputs are an alignment file which must contain
alignments for all possible queries matching the target sequence. This file
must be in [Stockholm](https://sonnhammer.sbc.su.se/Stockholm.html) format
with several custom metadata fields. The other required input is a set of
substitution matrices. This file uses a custom, but extremely simple format.
	
#### Alignment File Format

Alignments for all possible queries matching the target sequence should be
contained in a single file in Stockholm format.

There are several special metadata fields that must exist for each alignment in
this file. See the example below. An explanation is indented to the right of
each field with additional detail as noted.

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
#=GF GI -25                       gap init
#=GF GE -5                        gap extension
```

  * (1) query sequence names must be in the format 'name#family/class'
  * (2) valid values: 'q' if the alignment is on the reverse
    strand and the reversed sequence is the query; 't' if the alignment
    is on the reverse strand and the reversed sequence is the target;
    '-1' if the alignment is on the positive strand
  * (3) the flanking region of the unaligned query sequence
  * (4) the name of the substitution matrix file used to create alignment

#### Converting Alignments

PolyA can convert a Cross Match or Repeat Masker alignment file to the
particular version of Stockholm format it requires. To do this, use either the
`--cm-to-stockholm` or `--rm-to-stockholm` options, repspectively, passing the
path to the file to be converted.

The script will produce two files in the same directory as the input, one with a
`.sto` extension and the other with a `.matrix` extension. These can be passed
to the PolyA command line tool as `ALIGNMENTS` and `MATRICES`, respectively (see
`--help`).

Example:

```
python -m polyA --cm-to-stockholm my_alignments.cm
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

#### Sequence File (Optional)

A FASTA file of the target sequence is needed when using ULTRA.
The target sequence must be the same genomic region that was used 
to get the cross_match alignment file. This file must follow the
format of  
```
>chrom:start-end
target_sequence  
```
as shown in the example.

Sequence file example format:
```
>chr1:152302175-152325203
AATAGTTTATTTTTAATTTAGATGCAGCTTACTATAATATTAATTATGTCCAAGATGATT
TTTTGAATACAGAATACTAGAATTCCAATAGAAGGATAATAGAGAAAGATGTGCTAGCCC
...
```

### Output Formats

```
start   stop    ID  name
----------------------------------------
11990879    11991268    eaa042dd09f944f68dba2fd4727c64e2    LTR40a#LTR/ERVL
11991272    11991444    fb5ef5e0e2ca4e05837ddc34ca7ef9e4    MSTA1#LTR/ERVL-MaLR
11991445    11991562    bdfc4039b7d947d0b25bf1115cc282ed    AluJr4#SINE/Alu
11991563    11991573    4871d91441a146209b98f645feae68c8    FLAM_C#SINE/Alu
11991574    11991818    fb5ef5e0e2ca4e05837ddc34ca7ef9e4    MSTB1#LTR/ERVL-MaLR
11991819    11991875    eaa042dd09f944f68dba2fd4727c64e2    LTR40a#LTR/ERVL

* Matching IDnums correspond to partial sequences that originate from 
the same ancestral sequence.
```

#### Confidence only output file format

Computes confidence of a single input alignment region. Does not perform 
annotation or adjudication, simply outputs the confidence of all competing 
queries given in the input.

```
query_label         confidence
LTR40a#LTR/ERVL     0.875
LTR40b#LTR/ERVL     0.052
LTR40c#LTR/ERVL     0.001
...
```

### Extensions

#### Visualizing annotations using SODA

The command line option --soda will output the annotation data to a json file 
(output.0.viz) that can be used for visualization in SODA (linked below).
The json file can be submitted on the browser to view the TE annotations from PolyA
as well as the annotations from the UCSC Genome Browser for the same region of the
human genome (hg38). The PolyA visualization can display the confidence values for all competing 
annotations of a selected region as well as their corresponding sequence alignments.

https://sodaviz.cs.umt.edu/polya-soda.html

#### Prior Counts Files

Default confidence calculations assume a uniform distribution over all
competing queries. In the case of non uniform priors, the command line option --prior-counts prior_counts.txt includes prior
genome counts in confidence calculations (see paper for more details). 

https://www.biorxiv.org/content/10.1101/2021.02.13.430877v1

Prior counts file example format:
```
subfamily   genome_count
AluYk2      6855
LTR38	    255
L1PA7_5end  13261
...
```

#### Using ULTRA

The optional use of ULTRA allows polyA to include tandem repeats (TRs) in the competing annotations
of the target sequence. Doing so removes the dependency on pre-masking TRs prior to annotation, allows 
TRs to outcompete potentially weak fragmentary family annotation, and allows a family annotation
to outcompete a TR.
The command line option --sequences seq.fasta (with --ultra-path if necessary) will
run ULTRA with polyA or --ultra-data ultra_data.txt can be used if ULTRA was ran on seq.fasta prior.


### External software dependencies

#### Easel 
PolyA internally adjusts alignment scores to account for the scale multiplier
implicit to each score matrix (the so-called `lambda` value).
Unless the lambda value is given in the matrix file, PolyA must compute it; this
is done using the `esl_scorematrix` utility found in the 
[Easel](https://github.com/EddyRivasLab/easel/) software package. 

To prepare that tool, follow these steps:
```
    # get source code
    % cd $HOME/git  # or wherever you like to place repositories
    % git clone https://github.com/EddyRivasLab/easel
    % cd easel
    % autoconf
    % ./configure
    % make
    % gcc -g -Wall -I. -L. -o esl_scorematrix -DeslSCOREMATRIX_EXAMPLE esl_scorematrix.c -leasel -lm


    # add the easel directory to your path, e.g.
    % export PATH=$HOME/git/easel:$PATH    
    # you may wish to add this to your path permanently, e.g. 
    #    https://opensource.com/article/17/6/set-path-linux

    # alternatively, you can set the --easel-path argument to $HOME/git/easel
```

#### ULTRA

PolyA can optionally include tandem repeat annotations from ULTRA as competitors in
the adjudication process. Options are:

(i) If you have an ultra annotation output on hand for the sequence being annotated, reference
the ultra output file with the --ultra-data flag, e.g. 
```
    % polyA --ultra-data ultra.file ALIGNMENTS MATRICES
```

(ii) If you wish to run [ULTRA](https://github.com/TravisWheelerLab/ultra) on the fly, you'll 
need the software:
```   
    % cd $HOME/git  # or wherever you like to place repositories
    % git clone https://github.com/TravisWheelerLab/ULTRA ultra
    % cd ultra
    % cmake .
    % make
```
then you'll reference the ultra binary in the PolyA command:
```
% polyA --ultra-path $HOME/git/ultra/ultra --sequences seq.fasta  ALIGNMENTS MATRICES
```

## Development

This project uses [Pipenv](https://pipenv.pypa.io/en/latest/), which can be
installed through Homebrew for Mac users. It must be installed before the
Makefile targets, or the other commands listed in this document will work.

In order to run a command which relies on the project virtual environment, such
as `python foo.py`, it is necessary to either run `pipenv shell` first, which
will put you into a shell that has the correct Python in its `PATH`, or prefix
the command with `pipenv run` (e.g. `pipenv run python foo.py`).

To build a PyPI package, use `make build-package`. To publish, use
`make publish-package`. We use [Flit](https://flit.readthedocs.io/en/latest/) to
handle building and publishing a package, so it is also possible to invoke Flit
directly to do anything else it is capable of.

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

### Docker Images

There is a `Dockerfile` in the repo root. The image it describes is
used for running tests in CI and can be used locally for convenience.

When dependencies change the image must be rebuilt and the new
version pushed to Docker Hub. This can be done with `make container`
if you have the correct permissions. Otherwise, ask a maintainer
to do it for you.

There is also a file called `Dockerfile_run`. This is a runner image that is
available for users who prefer to run PolyA through Docker (which is often more
convenient). The runner image can be built with `tool/build-runner-image.sh`
and, assuming sufficient permissions, pushed to Docker Hub with
`tool/push-runner-image.sh`.

The scripts assume you are building the latest version of the image and so set
the `latest` tag. If this is not the case, it is fairly simple to issue the
correct Docker commands manually.
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

### Documentation

We use Sphinx for project documentation. Run `make docs` to update the
documentation sources and build the HTML. Use `make docs-serve` to serve the
documentation locally.

### Release Process

Generating a new release is a four step process. The first step MUST be
completed first, but the rest may be done in any order.

First, update the value of `VERSION` in `polyA/_version.py`, incrementing the
various portions of the version number accordingly. We would like to follow
[Semantic Versioning](https://semver.org), at least in general. Commit your
changes to the `master` branch.

Second, [create](https://github.com/TravisWheelerLab/polyA/releases/new) a
"release" on GitHub. The tag should consist of a "v", followed immediately by
the version number you chose above. For example: `v1.2.3`.

Third, run `make build-package`, followed by `make publish-package` assuming
there aren't any build errors.

Fourth, and finally, create and push a new runner image to Docker Hub. Build the
image by running `./tool/build-runner-image.sh`, and publish it with
`./tool/push-runner-image.sh`. The version tag will be set automatically and the
`latest` tag will also be updated.

## License

BSD license. See `LICENSE`.

## Authors

[Wheeler Lab](http://wheelerlab.org) at the University of Montana.

  - Kaitlin Carey
  - Audrey Shingleton
  - George Lesica
  - Jack Roddy
  - Travis Wheeler
