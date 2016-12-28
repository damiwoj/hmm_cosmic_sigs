#!/bin/bash

taxid=$1

SISTdir=~/supercoiling/Software/SIST

#usage: ./master.pl -f <sequence file> -a <algorithm_type> (choose algorithm_type: M, Z, C, or A) [options]
#
#Input requires a sequence file and algorithm type.
#Script analyzes superhelically induced structural transition probabilities for each base pair.
#Algorithm types: M (melting), Z (Z-DNA), C (cruciforms), and A (competition between all three).
#Sequence file will be converted to the format required by the algorithm.
#Code in directory trans_three/ will handle M, Z, and C algorithm types.
#Code in directory trans_compete/ will handle A algorithm type.
#For algorithm type options -a C and -a A user will need an Inverted Repeat Finder (IRF) executable compatible with user's operating system.
#IRF download page: http://tandem.bu.edu/irf/irf.download.html.
#Selected output will be printed to the screen.
#
#Options:
#-f       Required: specify sequence file 
#-a       Required: specify algorithm type: M, Z, C, or A
#-t       Optional: set temperature (default 310)
#-s       Optional: set superhelical density (default -0.06)
#-i       Optional: set ionic strength (default 0.01)
#-th      Optional: set energy threshold (default 12), -t 10 is recommended for -a A (competition algorithm)
#-c       Optional: flag to set molecular type to circular (default linear)
#-n       Optional: flag to set melting energetics to nearest neighbor (default copolymeric)
#-b       Optional: flag to print base pair for each position (default null)
#-p       Optional: flag to print algorithm parameters (default null)
#-r       Optional: flag to print ensemble average results (default null)


function do_it {
    file1=$1
    file2=$2
    curdir=`pwd`
    cd ${SISTdir}
    ./master.pl -f ${curdir}/${file1}.seq -a A -s -0.05 -c > ${curdir}/${file1}.txt
    cd ${curdir}
    python generate_next_cycle.py ${file1}.seq ${file1}.txt 0.01 1> ${file2}.seq 2> ${file1}.sibz
}

do_it 00 01
do_it 01 02
do_it 02 03
do_it 03 04
do_it 04 05
do_it 05 06
