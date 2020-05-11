#!/usr/bin/env bash

# perl AdjudicateRegions.pl --lambda .1227 ex3_hg38_nesting.align.format 25p41g_edited.matrix | tee perl-output.txt
# echo ""
# python3 AdjudicateRegions.py --matrixPos --lambda .1227 --matrixPos Seqs_anotherExample.align 25p41g_edited.matrix | tee python-output.txt

for f in *.align ;
	do echo $f;
		python3 AdjudicateRegions.py --lambda .1227 $f 25p41g_edited.matrix;
		echo "-------------------------------------------------------------\n";
	done

for f in *.align.format ;
	do echo $f;
		python3 AdjudicateRegions.py --lambda .1227 $f 25p41g_edited.matrix;
		echo "-------------------------------------------------------------\n";
	done


# if [ `uname -s` == Darwin ]; then
#   DATECMD=gdate
# else
#   DATECMD=date
# fi
# 
# ## ---------
# ## UTILITIES
# ## ---------
# 
# function compare-elapsed() {
#   local base_file="$1"
#   local pl_elapsed=`cat "$base_file.pl.elapsed"`
#   local py_elapsed=`cat "$base_file.py.elapsed"`
#   echo -e "Perl: $pl_elapsed, Python: $py_elapsed - $base_file"
# }
# 
# function run-script() {
#   local alignment="$1"
#   local interpreter="$2"
#   local extension="$3"
#   local start=`$DATECMD +'%s%N'`
#   "$interpreter" "AdjudicateRegions.$extension" \
#       --lambda 0.1227 \
#       "$alignment" \
#       25p41g_edited.matrix \
#       > "$alignment.$extension.output" \
#       2>&1
#   local finish=`$DATECMD +'%s%N'`
#   local delta=`expr \( $finish / 1000000 \) - \( $start / 1000000 \)`
#   echo "$delta" > "$alignment.$extension.elapsed"
# }
# 
# function run-perl() {
#   run-script "$1" "perl" "pl"
# }
# 
# function run-python() {
#   run-script "$1" "python" "py"
# }
# 
# ## -----
# ## TESTS
# ## -----
# 
# # artificial test files most often used
# for alignment in *.align; do
#   run-perl "$alignment"
#   run-python "$alignment"
#   compare-elapsed "$alignment"
#   # assert-same "$alignment"
# done
# 
# # real test files from ISB, not as frequently used but need to use later 
# for alignment in *.align.format; do
#   run-perl "$alignment"
#   run-python "$alignment"
#   compare-elapsed "$alignment"
#   # assert-same "$alignment"
# done
