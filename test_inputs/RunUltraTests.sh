#!/usr/bin/env bash

# check for valid input
if [ "$#" -ne 2 ]
then
  echo "$0: takes two arguments: <test__file_directory> <path_to_ultra>"
  exit 1
fi

test_file_dir=$1
ultra_path=$2
test_files=$(find $test_file_dir -type f -name "*.fa")

# run AdjudicateRegions
for file in $test_files
  do
    echo $file
    # alignment files will correspond with fasta file names from CM
  	align_file="${file}.cm.align"
  	echo "$align_file"
  	poetry run python3 AdjudicateRegions.py --seqpos --lambda .1227 --ultraPath $ultra_path --seqFile $file $align_file 25p41g_edited.matrix
    echo "-------------------------------------------------------------"
done
	