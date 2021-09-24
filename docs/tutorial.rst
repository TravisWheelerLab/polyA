Tutorial
========

This tutorial covers running PolyA against region chr2:128348374-128348687.

Input Data
----------

#. Get the FASTA file for the target region (chr2_128348379_128354757.fa)

#. Obtain the cross_match alignment file for this region (chr2_128348379_128354757.fa.cm)

#. In the polyA directory, convert the cross_match file to stockholm format and get the required substitution matrices
    * % polyA --cm-to-stockholm chr2_128348379_128354757.fa.cm
    * This will output the following files needed to run polyA:
         * chr2_128348379_128354757.fa.cm.sto
         * chr2_128348379_128354757.fa.cm.matrix

The files used here can found in the polyA tutorial directory.

Generating Output
-----------------

Run polyA:
    % polyA --sequence-position chr2_128348379_128354757.fa.cm.sto chr2_128348379_128354757.fa.cm.matrix

If you wish to include tandem repeats, you can run ULTRA on the fly:
    % polyA --ultra-path /path/to/ultra --sequences chr2_128348379_128354757.fa chr2_128348379_128354757.fa.cm.sto chr2_128348379_128354757.fa.cm.matrix

or you can run ULTRA in advance:
    % ./ultra -ss chr2_128348379_128354757.fa > chr2_128348379_128354757.fa.ultra
and then run polyA with the following command:
    % polyA --sequence-position --ultra-data chr2_128348379_128354757.fa.ultra chr2_128348379_128354757.fa.cm.sto chr2_128348379_128354757.fa.cm.matrix
