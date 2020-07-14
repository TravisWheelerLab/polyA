#!/usr/bin/env bash

if [ `uname -s` == Darwin ]; then
  DATECMD=gdate
else
  DATECMD=date
fi

python3 AdjudicateRegions.py --lambda .1227 ex3_hg38_nesting.align.format 25p41g_edited.matrix

# for f in *.align ;
# 	do echo $f;
# 		python3 AdjudicateRegions.py --lambda .1227 $f 25p41g_edited.matrix;
# 		echo "-------------------------------------------------------------\n";
# 	done

# for f in *.align.format ;
# 	do echo $f;
# 		python3 AdjudicateRegions.py --lambda .1227 $f 25p41g_edited.matrix;
# 		echo "-------------------------------------------------------------\n";
# 	done
