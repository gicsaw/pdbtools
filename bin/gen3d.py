#!/usr/bin/env python
import subprocess
import pdbtools.ligand_tools as ligand_tools


def main():

    import argparse
    title_line = 'generate 3d conformation from SMILES'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_smi', required=True,
                        help='input SMILES')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output ligand pdb file')

    args = parser.parse_args()
    smi = args.input_smi
    output_file = args.output_file

    e = ligand_tools.gen_3d(smi, output_file)
    if e is not None:
        print(e)

if __name__ == "__main__":
    main()
