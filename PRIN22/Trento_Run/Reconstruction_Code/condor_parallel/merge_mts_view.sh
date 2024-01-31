#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Script for improving the micro-tracks after linking. Usage is: '
    echo ' '
    echo 'source merge_mts_view.sh viewID'
    echo ' '
    return 0
fi

python3 merge_mts_view.py $1
#echo $1
