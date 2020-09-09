# Perl to Python Tests

The script `RunTests.sh` contains tests to compare the output of the Python
rewrite with the output of the original Perl script. Since the Python script
has evolved, these are probably not useful, so these comparisons are
commented out. The script also contains tests to compare the runtime of the
two versions of the code. This is what runs by default.

The following procedure will run the tests:

```
cd polyA
pipenv install
brew install coreutils # if on a Mac
cd test_inputs
./RunTests.sh
```

Sample output recorded on a MacBook Air, times are in ms:

```
Perl: 1092, Python: 36 - L1M5_orf2.10.align
Perl: 1457, Python: 37 - Seqs_AluAlu.align
Perl: 94, Python: 34 - Seqs_AluAlu_subset.align
Perl: 2298, Python: 38 - Seqs_NoNesting.align
Perl: 4998, Python: 40 - Seqs_anotherExample.align
Perl: 4321, Python: 40 - Seqs_doubleNested.align
Perl: 231, Python: 34 - Seqs_doubleNested_subset.align
Perl: 5412, Python: 41 - Seqs_doubleNesting2.align
Perl: 191, Python: 33 - Seqs_gap.align
Perl: 4603, Python: 40 - Seqs_noStitching.align
Perl: 4398, Python: 39 - seqs_fullAlu.align
Perl: 332, Python: 35 - seqs_fullAlu_subset2.align
Perl: 9261, Python: 40 - ex10_hg38_alugedon.align.format
Perl: 6251, Python: 38 - ex1_hg38_nesting.align.format
Perl: 18596, Python: 40 - ex2_hg38_nesting.align.format
Perl: 1369, Python: 35 - ex3_hg38_nesting.align.format
Perl: 115, Python: 35 - ex4_hg38_badoverlap.align.format
Perl: 26987, Python: 42 - ex5_hg38_composite.align.format
Perl: 4871, Python: 44 - ex6_hg38_ambiguous.align.format
Perl: 3175, Python: 36 - ex7_hg38_ambiguous.align.format
Perl: 11502, Python: 40 - ex8_hg38_alugedon.align.format
Perl: 5554, Python: 39 - ex9_hg38_alugedon.align.format
```
