#!/bin/bash
list_file=$1
mkdir pdb
line_pdb="wget"
line_fasta="wget"
while read line
do
    array=($line)
    #protein=${array[0]}
    code=${array[0]}
    echo $code
    line_pdb+=" https://files.rcsb.org/download/$code.pdb"
#    wget https://www.rcsb.org/fasta/entry/$code -O pdb/$code.fasta
    line_fasta+=" https://www.rcsb.org/fasta/entry/$code"

done < $list_file

$line_pdb -P pdb
$line_fasta -P pdb
while read line
do
    array=($line)
    #protein=${array[0]}
    code=${array[0]}
    mv pdb/$code  pdb/$code.fasta

done < $list_file

