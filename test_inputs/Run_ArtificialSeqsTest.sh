  #!/usr/bin/env bash

if [ `uname -s` == Darwin ]; then
  DATECMD=gdate
else
  DATECMD=date
fi

## ---------
## UTILITIES
## ---------

# for f in /Users/kaitlincarey/Desktop/NOTEBOOK/2020/Adjudication/Aug25-TestArtificalSeqs/TEST_ARTIFICALSEQS/*/*.align;
# for f in /Users/kaitlincarey/Desktop/NOTEBOOK/2020/Adjudication/Aug25-TestArtificalSeqs/TEST_ALU_SEQS_NORECOMB/*.align;
 for f in /Users/kaitlincarey/Desktop/NOTEBOOK/2020/Adjudication/Aug25-TestArtificalSeqs/TEST_ARTIFICIALSEQS_NESTING/*.align;
 	do echo $f;
 		poetry run python3 AdjudicateRegions.py --seqpos --lambda .1227 $f 25p41g_edited.matrix;
 		echo "-------------------------------------------------------------\n";
 	done

