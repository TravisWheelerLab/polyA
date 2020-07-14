  #!/usr/bin/env bash

if [ `uname -s` == Darwin ]; then
  DATECMD=gdate
else
  DATECMD=date
fi

## ---------
## UTILITIES
## ---------

 for f in *.align.format ;
 	do echo $f;
 		poetry run python3 AdjudicateRegions.py --seqpos --lambda .1227 $f 25p41g_edited.matrix;
 		echo "-------------------------------------------------------------\n";
 	done

