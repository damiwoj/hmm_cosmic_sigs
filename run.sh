#!/bin/bash

taxid=$1



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
