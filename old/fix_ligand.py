#!/usr/bin/env python
from pifinder.pdbtools import pdbtools


def main():

    import argparse
    title_line = 'Fixer for ligand pdb'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', type=str, required=True,
                        help='input ligand pdb file')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='o.pdb',  help='output ligand pdb file')
    parser.add_argument('-n', '--neutralize', action='store_true',
                        required=False, help='neutralize mol ')
    parser.add_argument('-p', '--pH', type=float, required=False, default=None,
                        help='protonation state')
    parser.add_argument('--is_peptide', action='store_true',
                        required=False, help='is peptide?')
    parser.add_argument('--bond_rebuilding', action='store_true',
                        required=False, help='bond_order_rebuilding?')
    parser.add_argument('--fix_atom_idx', action='store_true',
                        required=False, help='fix_atom_idx')

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    neutralize = args.neutralize
    pH = args.pH
    is_peptide = args.is_peptide

    bond_rebuilding = args.bond_rebuilding
    is_fix_atom_idx = args.fix_atom_idx
    out_format = output_file.split('.')[-1]
    if out_format == 'pdbqt':
        tmp_file = '.'.join(output_file.split('.')[:-1]) + '.pdb'
    else:
        tmp_file = output_file
    pdbtools.fix_ligand(input_file, output_file, neutralize=neutralize,
                        pH=pH, is_peptide=is_peptide,
                        bond_rebuilding=bond_rebuilding,
                        is_fix_atom_idx=is_fix_atom_idx)
    if out_format == 'pdbqt':
        pdbtools.ligand_to_pdbqt(tmp_file, output_file)


if __name__ == "__main__":
    main()
