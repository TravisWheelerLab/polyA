#!/bin/sh

# artificial test files most often used
#for f in *.align ;
#do echo $f;
    #perl AdjudicateRegions.pl --lambda .1227 $f 25p41g_edited.matrix;
    #python3 AdjudicateRegions.py --lambda .1227 $f 25p41g_edited.matrix;
    #echo "-------------------------------------------------------------\n";
#done

perl AdjudicateRegions.pl --lambda .1227 seqs_fullAlu_subset2.align 25p41g_edited.matrix | tee perl-output.txt
python3 AdjudicateRegions.py --lambda .1227 seqs_fullAlu_subset2.align 25p41g_edited.matrix | tee python-output.txt

# real test files from ISB, not as frequently used but need to use later 
#for f in *.align.format ;
#do echo $f;
#perl AdjudicateRegions.pl --lambda .1227 $f 25p41g_edited.matrix;
#echo "-------------------------------------------------------------\n";
#done