#!/usr/bin/env python
from pdbtools.pdbtools import pdbtools


def main():

    import argparse
    title_line = 'Fixer for protein pdb which is converted from pdbqt'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input protein pdb file')
    parser.add_argument('-r', '--ref_file', type=str, required=True,
                        help='reference protein pdb file')
    parser.add_argument('-t', '--tmp_dir', type=str, required=False,
                        default='.', help='temp directory, defalut=.')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdb', help='output protein pdb file')

    args = parser.parse_args()
    input_file = args.input_file
    ref_file = args.ref_file
    tmp_dir = args.tmp_dir

    output_file = args.output_file

    result = pdbtools.read_pdb_protein(input_file, remain_remark=True)
    (protein_chain_residues, pdb_info_lines, protein_dict,
     ligand_dict, conect_dict) = result

    result = pdbtools.align_3d(input_file, protein_dict,
                               ligand_dict, ref_file, tmp_dir)
    (protein_dict_new, ligand_dict_new) = result
    molecule_dict = dict()
    molecule_dict.update(protein_dict_new)
    molecule_dict.update(ligand_dict_new)

    model_dict = dict()
    model_dict[0] = (pdb_info_lines, molecule_dict, conect_dict)
    pdbtools.write_model_pdb(model_dict, output_file)


if __name__ == "__main__":
    main()
