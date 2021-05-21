#!/bin/bash

list_file=$1
ref_file=$2
while read line
do
    array=($line)
    #protein=${array[0]}
    chain_file=${array[0]}
#    code=${array[0]}
    echo $chain_file
    split_ligand.py -i $chain_file -r $ref_file -e HOH -d chain/

done < $list_file
#ls ????A_???.pdb > list_ligand.txt

