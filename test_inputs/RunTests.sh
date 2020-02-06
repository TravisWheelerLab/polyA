#!/bin/sh

# artificial test files most often used
for f in *.align ;
do echo $f;
perl AdjudicateRegions.pl --lambda .1227 $f 25p41g_edited.matrix;
echo "-------------------------------------------------------------\n";
done


# real test files from ISB, not as frequently used but need to use later 
for f in *.align.format ;
do echo $f;
perl AdjudicateRegions.pl --lambda .1227 $f 25p41g_edited.matrix;
echo "-------------------------------------------------------------\n";
done
