#!/bin/bash
# run 201, 1K and All, info overview analysis on the scafSeq file of argument $1
# results are append to the file prefixed by arg $2
# this file name is $2_scafContigInfo.txt
# you usually use asm and a number, for example asm28
#
# an example call:
#      scafSeqContigInfo.sh *scafSeq asm28
if [ $# -ne 2 ]; then
	echo
    echo "scafSeqContigInfo.sh does runs of the script to compute scaffold info using several scaffold size cutoffs."
	echo "It takes the scaffold fasta file and an assembly prefix for the output file as arguments, e.g."
	echo
	echo "      scafSeqContigInfo.sh *scafSeq asm28"
	echo "         # this example would output to a file named asm28_scafContigInfo.txt"
	echo
else
	out=${2}_scafContigInfo.txt
	echo $1 > $out
	scafSeqContigInfo.py $1 -x -S 300 >> $out
	echo >> $out
	scafSeqContigInfo.py $1 -nc -x -s 201 >> $out
	echo >> $out
	scafSeqContigInfo.py $1 -nc -x -s 1K >> $out
	echo >> $out
	scafSeqContigInfo.py $1 -nc -x -All  >> $out
fi
