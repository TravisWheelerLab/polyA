![PolyA CI](https://github.com/TravisWheelerLab/polyA/workflows/PolyA%20CI/badge.svg)

# AAAAAAAAAAAAAAAA (PolyA)

> **A**utomatically **A**djudicate **A**ny **A**nd **A**ll **A**rbitrary
> **A**nnotations, **A**stutely > **A**djoin **A**butting **A**lignments,
> **A**nd **A**lso **A**mputate **A**nything **A**miss.

## About

**TODO:** Kaitlin

A common annotation process compares an unannotated sequence to a collection of known sequences.
Sometimes, more than one of the queries shows a significant match to the target. Current 
annotation processes choose the 'true' match based on highest alignment score. In databases
with highly similar sequences, more than one query may match with high alignment scores, this
method is no longer reliable because it is possible that either query is the true sequence, and 
it falsely implies certainty in the true match. PolyA produces confidence estimates for 
competing annotations to eliminate implying false certainty in the annotations. In addition, 
polyA can identify instances of gene conversion and homologous recombination, and find the 
exact genomic location of the switch between sequences. It can also clarify ambiguous boundaries
between neighboring elements due to homologous over extension. And it can find sequences nesting
while identifying the inserted sequence along with the original sequence that was inserted into.


## The Algorithm

**TODO:** Kaitlin - we can also link to whatever paper(s) are published

Using our confidence score analysis enables a jumping HMM approach, allowing transitions between
competing queries identifies switches between queries. These switches are a result of gene
conversion, homologous recombination, neighboring elements and sequence nesting.

Our algorithm has 4 parts: 

1. Confidence calculations help choose between competing annotations
	using an application of the bayes theorem on the probabilities underlying alignments scores
	we compute posterior probabilities (confidence) of the likelihood a query is the true
	source for the genomic region.
2A. Segmented adjacent confidence finds non-exact breakpoint locations
	-takes alignments of all competing queries
	-splits alignments up into contiguous segments
	-compute confidence scores for segments 
	-dynamic programming pass through matrix holding confidence values gives most probable
	path through alignment
	-allows switches between queries
	-identifies a query for each segment
2B. Segmented overlapping confidence finds exact breakpoints of gene conversion, homologous 
recombination and boundary detection
	-takes alignments of all competing queries
	-splits alignments up into adjacent segments starting at every nucleotide
	-compute average confidence values for each nucleotide based on all segments overlapping
	that position
	-new dynamic programming matrix with confidence values for each nucleotide instead
	of each segment
	-dynamic programming pass identifies a query for every nucleotide
3. Graph Algorithm finds nested sequences
	-all sequences labeled in dynamic programming become nodes
	-find alternative paths through that graph based on confidence of the sink subfamily
	in the source node
	-splice out all sequences with no alternative paths (these are the inserted sequences)
	-splice corresponding positions out of dynamic programming matrix
	-another dynamic programming pass on the updated matrix stitches original sequence back 
	together
	
For a more detailed description view the [poster](/supporting_materials/AlgorithmPoster.pdf)

## Using

**TODO:** Kaitlin, George

### Input File Format
	
#### Alignment Files

  - RM of CM alignment files 
  - Links to RM and CM alignment file formats
  - Add hmmer eventually?

Score matrix files example format:

```
  A   R   G   C   Y   T   K   M   S   W   N   H   V   X
  8   1  -6 -13 -14 -15 -11  -2 -10  -3  -1  -1  -1 -27
  2   2   1 -13 -13 -14  -6  -5  -5  -5  -1  -1  -1 -27
 -2   3  10 -13 -13 -13  -1  -8  -1  -8  -1  -1  -1 -27
-13 -13 -13  10   3  -2  -8  -1  -1  -8  -1  -1  -1 -27
-14 -13 -13   1   2   2  -5  -6  -5  -5  -1  -1  -1 -27
-15 -14 -13  -6   1   8  -2 -11 -10  -3  -1  -1  -1 -27
 -9  -5  -1  -9  -6  -2  -1  -9  -5  -5  -1  -1  -1 -27
 -2  -6  -9  -1  -5  -9  -9  -1  -5  -5  -1  -1  -1 -27
 -8  -4  -1  -1  -4  -8  -4  -4  -1  -8  -1  -1  -1 -27
 -3  -6 -10 -10  -6  -3  -6  -6 -10  -3  -1  -1  -1 -27
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -27
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -27
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -27
-27 -27 -27 -27 -27 -27 -27 -27 -27 -27 -27  -1  -1 -27
```

### Output file format

#### Genomic Location

```
start stop	IDnum	query
0	362		1111	L1PREC2_3end
363	567	2345	AluJr
568	833	3579	AluYb8
834	964	1245	AluJr
965	980	6047	L1MA4A_3end
981	1497	1111	L1PREC2_3end
```

**FIXME** - switch these postions to genomic locations not matrix pos 

Matching IDnums correspond to sequences involved in nesting that have been
stitched back to the original sequence.


### Additional software

esl_scorematrix as a part of the esl package in the hammer software suite
	-will compute lambda for the input score matrix
	-not needed if including lambda as a command line argument

### Using at the command line

```
usage: $0 alignFile matrixFile
ARGUMENTS
	--gapInit [-25]
	--getExt [-5]
	--lambda [will calc from matrix if not included - need esl_scorematrix installed]
	--segmentsize [30]
	--changeprob [1e-45]
	
OPTIONS
	--help - display help message
	--printmatrices - output all dynamic programming matrices
	--matrixpos - prints subfam changes in matrix position instead of genomic position
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

### Perl-to-Python Tests

```
pipenv install # in the repo root
pipenv shell
```

Once you are in the virtual environment shell, you can run
`PYTHONPATH=../ python ./RunTests.sh` in the `test_inputs` directory.
This is also possible using `make check-slow` and happens in CI.

## License

MIT license. See `LICENSE`.

## Authors

[Wheeler Lab](http://wheelerlab.org) at the University of Montana.

  - Kaitlin Carey
  - Audrey Shingleton
  - George Lesica
  - Travis Wheeler
