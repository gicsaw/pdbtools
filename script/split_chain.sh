#!/bin/bash

list_file=$1
mkdir chain
while read line
do
    array=($line)
    #protein=${array[0]}
    code=${array[0]}
    echo $code
    split_chain.py -i pdb/$code.pdb -d chain/

done < $list_file
#ls ?????.pdb 

