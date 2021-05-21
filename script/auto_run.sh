#!/bin/bash
### find target protein at uniprot https://www.uniprot.org
### See "Structure" in uniprot website
### Collect PDB ID by referring to position information and save it in list.txt 
./dw.sh list.txt
./split_chain.sh list.txt
ls chain/?????.pdb > list_chain.txt
### choose reference pdb (Ex: 2P2IA) for alignment 
./split_ligand.sh list_chain.txt chain/3HMMA.pdb
ls chain/?????_???.pdb > list_ligand.txt
### if you want only chain A, ls chain/????A_???.pdb > list_ligand_A.txt
### or grep A_ list_ligand.txt > list_ligand_A.txt
### choose reference ligand (Ex: 2P2IA_608) to fine pocket ligands
./dist.sh list_ligand.txt chain/3HMMA_855.pdb > list_ligand_new.txt
cp list_ligand_new.txt list_ligand_select.txt
./pymol.sh list_ligand_new.txt
### run_output script and check ligands
### fix list_ligand_select.txt  manualy 
### vi list_ligand_select.txt
./cp_pdb.sh list_ligand_select.txt > list_final.txt
### see select directory
./dw_lig_ref.sh list_final.txt

./fix.sh list_final.txt

./auto_box.sh list_final.txt

