#!/bin/bash

# This bach script replaces the atoms section inside of the conf_lammps.in file written by icnf.exe

if [ $# -ne 3 ]; then
    echo "Usage: $0 <configuration file from icnf.exe> <new atom coordinates> <output configuration file>"
    exit 1
fi
CONFFILE=$1
NEWCOORDS=$2
OUTPUTFILE=$3

START=`grep -n Atoms ${CONFFILE} | cut -f 1 -d ':'`
END=`grep -n Bonds ${CONFFILE} | cut -f 1 -d ':'`
let START=$START-1
let END=$END-1

cmd="sed -n '1,$START p'<${CONFFILE} > tmp.start"
eval $cmd 
cmd="sed -n '$END,$ p'<${CONFFILE} > tmp.end"
eval $cmd 

cat tmp.start ${NEWCOORDS} tmp.end > ${OUTPUTFILE}
