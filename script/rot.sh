#!/bin/bash

list_file=$1
ref_file=$2
while read line
do
    array=($line)
    #protein=${array[0]}
    code=${array[0]}
    echo $code
    align_3d.py -i chain/$code.pdb -r chain/$ref_file -o chain/"$code"_rotate.pdb

done < $list_file

