About
=====

Preprint on bioRxiv `<https://www.biorxiv.org/content/10.1101/2021.02.13.430877v1/>`_

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