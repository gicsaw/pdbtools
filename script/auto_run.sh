#!/bin/bash
### find target protein at uniprot https://www.uniprot.org
### See "Structure" in uniprot website
### Collect PDB ID by referring to position information and save it in list.txt 

uniprot_dict.py pdb EGFR_HUMAN
filter_pdb_list.py --resolution 10 -r 650- P00533_pdb_list.txt > P00533_pdb_list_filter.txt
dw_pdb.py P00533_pdb_list_filter.txt pdb
find_mutation_pdb.py P00533_pdb_list_filter.txt P00533.fasta pdb > list_mutation.txt

./split_chain.sh list_mutation.txt
python gen_chain_list.py -i list_mutation.txt -o list_chain.txt
### choose reference pdb (Ex: 2P2IA) for alignment 
./rot.sh list_chain.txt chain/3W2SA.pdb
./split_ligand.sh list_chain.txt
ls chain/?????_???.pdb > list_ligand.txt
### if you want only chain A, ls chain/????A_???.pdb > list_ligand_A.txt
### or grep A_ list_ligand.txt > list_ligand_A.txt
### choose reference ligand (Ex: 2P2IA_608) to fine pocket ligands
./dist.sh list_ligand.txt chain/2ITYA_IRE.pdb > list_ligand_new.txt
cp list_ligand_new.txt list_ligand_select.txt
./pymol.sh list_ligand_new.txt
### run_output script and check ligands
### fix list_ligand_select.txt  manualy 
### vi list_ligand_select.txt
./cp_pdb.sh list_ligand_select.txt > list_final.txt
./check_water.sh list_ligand_select.txt
ls select/*_HOH.pdb > water_file_list.txt
python consensus_water_coor.py list_ligand_select.txt > consensus_water_list.txt
#select water in consensus_water_list and remove other waters
python find_consensus_water.py select/1M17A_HOH.pdb consensus_water_list.txt

### see select directory

./dw_lig_ref.sh list_final.txt
./fix.sh list_final.txt
./auto_box.sh list_final.txt

# fix_protein.py -i select/"$code"_receptor.pdb -o fix/"$code"_receptor.pdb

