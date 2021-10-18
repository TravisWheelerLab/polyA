Running
=======

Run PolyA from the command line:

::

    polyA -h

or

::

    python -m polyA -h


Command line usage is also included below for convenience.

::

    usage: polyA [-h] [-v] [--chunk-size CHUNK_SIZE] [--trans-penalty TRANS_PENALTY]
              [--confidence] [--prior-counts FILE] [--shard-gap SHARD_GAP]
              [--sequences SEQS] [--ultra-data FILE] [--easel-path BIN]
              [--ultra-path BIN] [--log-file LOG] [--log-level LEVEL]
              [--matrix-position] [--output-path PATH] [--sequence-position] [--soda]
              ALIGNMENTS MATRICES

    positional arguments:
      ALIGNMENTS              alignments file in Stockholm format

      MATRICES                substitution matrices file in PolyA matrix format

    optional arguments:
      -h, --help              show this help message and exit
      -v, --version           show version and exit
      --chunk-size CHUNK_SIZE
                              size of the window in base pairs analyzed together
      --trans-penalty TRANS_PENALTY
                              penalty for changing annotations
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


Input Formats
-------------

PolyA accepts two required inputs and several optional inputs that affect
its behavior. The required inputs are an alignment file which must contain
alignments for all possible queries matching the target sequence. This file
must be in `Stockholm <https://sonnhammer.sbc.su.se/Stockholm.html>`_ format
with several custom metadata fields. The other required input is a set of
substitution matrices. This file uses a custom, but extremely simple format.

Alignment File
^^^^^^^^^^^^^^

Alignments for all possible queries matching the target sequence should be
contained in a single file in Stockholm format.

There are several special metadata fields that must exist for each alignment in
this file. See the example below. An explanation is indented to the right of
each field with additional detail as noted.

::

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

1. query sequence names must be in the format 'name#family/class'
2. valid values: 'q' if the alignment is on the reverse
   strand and the reversed sequence is the query; 't' if the alignment
   is on the reverse strand and the reversed sequence is the target;
   '-1' if the alignment is on the positive strand
3. the flanking region of the unaligned query sequence
4. the name of the substitution matrix file used to create alignment

Converting Alignments
^^^^^^^^^^^^^^^^^^^^^

PolyA can convert a Cross Match or Repeat Masker alignment file to the
particular version of Stockholm format it requires. To do this, use either the
``--cm-to-stockholm`` or ``--rm-to-stockholm`` options, respectively, passing the
path to the file to be converted.

The script will produce two files in the same directory as the input, one with a
``.sto`` extension and the other with a ``.matrix`` extension. These can be passed
to the PolyA command line tool as ``ALIGNMENTS`` and ``MATRICES``, respectively (see
``--help``).

Example:

::

    python -m polyA --cm-to-stockholm my_alignments.cm


Substitution Matrix File
^^^^^^^^^^^^^^^^^^^^^^^^

The substitution matrix file example format (can include ambiguity codes):

* this file must include all of the matrices specified in the
  ``#=GF MX`` field of the alignment file, with corresponding
  and matching matrix names
* if lambda is not included polyA will use esl_scorematrix to
  calculate it for all matrices

::

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


Sequence File
^^^^^^^^^^^^^

A FASTA file of the target sequence is needed when using ULTRA.
The target sequence must be the same genomic region that was used
to get the cross_match alignment file. This file must follow the
format of

::

    >chrom:start-end
    target_sequence

as shown in the example below.

::

    >chr1:152302175-152325203
    AATAGTTTATTTTTAATTTAGATGCAGCTTACTATAATATTAATTATGTCCAAGATGATT
    TTTTGAATACAGAATACTAGAATTCCAATAGAAGGATAATAGAGAAAGATGTGCTAGCCC
    ...


Output Formats
--------------

::

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


Confidence Output Format
^^^^^^^^^^^^^^^^^^^^^^^^

Computes confidence of a single input alignment region. Does not perform
annotation or adjudication, simply outputs the confidence of all competing
queries given in the input.

::

    query_label         confidence
    LTR40a#LTR/ERVL     0.875
    LTR40b#LTR/ERVL     0.052
    LTR40c#LTR/ERVL     0.001
    ...

Extensions
----------

Visualizing Annotations
^^^^^^^^^^^^^^^^^^^^^^^

The command line option ``--soda`` will output the annotation data to a json file
(output.0.viz) that can be used for visualization in SODA (linked below).
The json file can be submitted on the browser to view the TE annotations from PolyA
as well as the annotations from the UCSC Genome Browser for the same region of the
human genome (hg38). The PolyA visualization can display the confidence values for all competing
annotations of a selected region as well as their corresponding sequence alignments.

https://sodaviz.cs.umt.edu/polya-soda.html

Prior Counts Files
^^^^^^^^^^^^^^^^^^

Default confidence calculations assume a uniform distribution over all
competing queries. In the case of non uniform priors, the command line option --prior-counts prior_counts.txt includes prior
genome counts in confidence calculations (see paper for more details).

https://www.biorxiv.org/content/10.1101/2021.02.13.430877v1

Prior counts file example format:

::

    subfamily   genome_count
    AluYk2      6855
    LTR38	    255
    L1PA7_5end  13261
    ...

Using ULTRA
^^^^^^^^^^^

The optional use of ULTRA allows polyA to include tandem repeats (TRs) in the competing annotations
of the target sequence. Doing so removes the dependency on pre-masking TRs prior to annotation, allows
TRs to outcompete potentially weak fragmentary family annotation, and allows a family annotation
to outcompete a TR.
The command line option ``--sequences seq.fasta`` (with ``--ultra-path`` if necessary) will
run ULTRA with polyA or ``--ultra-data ultra_data.txt`` can be used if ULTRA was ran on seq.fasta prior.
