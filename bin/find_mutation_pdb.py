#!/usr/bin/env python
import sys
import os
from pdbtools import nwalign

usage = '''
mut_pdb.py list_file fasta_file pdb_dir
'''

Res31 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
         'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
         'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
         'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
         'ASX': 'N', 'GLX': 'Q', 'UNK': 'X', 'INI': 'K', 'AAR': 'R',
         'ACE': 'X', 'ACY': 'G', 'AEI': 'T', 'AGM': 'R', 'ASQ': 'D',
         'AYA': 'A', 'BHD': 'D', 'CAS': 'C', 'CAY': 'C', 'CEA': 'C',
         'CGU': 'E', 'CME': 'C', 'CMT': 'C', 'CSB': 'C', 'CSD': 'C',
         'CSE': 'C', 'CSO': 'C', 'CSP': 'C', 'CSS': 'C', 'CSW': 'C',
         'CSX': 'C', 'CXM': 'M', 'CYG': 'C', 'CYM': 'C', 'DOH': 'D',
         'EHP': 'F', 'FME': 'M', 'FTR': 'W', 'GL3': 'G', 'H2P': 'H',
         'HIC': 'H', 'HIP': 'H', 'HTR': 'W', 'HYP': 'P', 'KCX': 'K',
         'LLP': 'K', 'LLY': 'K', 'LYZ': 'K', 'M3L': 'K', 'MEN': 'N',
         'MGN': 'Q', 'MHO': 'M', 'MHS': 'H', 'MIS': 'S', 'MLY': 'K',
         'MLZ': 'K', 'MSE': 'M', 'NEP': 'H', 'NPH': 'C', 'OCS': 'C',
         'OCY': 'C', 'OMT': 'M', 'OPR': 'R', 'PAQ': 'Y', 'PCA': 'Q',
         'PHD': 'D', 'PRS': 'P', 'PTH': 'Y', 'PYX': 'C', 'SEP': 'S',
         'SMC': 'C', 'SME': 'M', 'SNC': 'C', 'SNN': 'D', 'SVA': 'S',
         'TPO': 'T', 'TPQ': 'Y', 'TRF': 'W', 'TRN': 'W', 'TRO': 'W',
         'TYI': 'Y', 'TYN': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y',
         'YOF': 'Y', 'FOR': 'X', '---': '-', 'PTR': 'Y', 'LCX': 'K',
         'SEC': 'D', 'MCL': 'K', 'LDH': 'K', 'CY0': 'C'}

Res20 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

Res13 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
         'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
         'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
         'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
         'X': 'UNK'}

Amino_type = "ARNDCQEGHILKMFPSTWYVX"
Ntype = 21
Amino_dict = dict()
for i in range(0, Ntype):
    amino_acid = Amino_type[i]
    Amino_dict[amino_acid] = i


def read_pdb(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()
    chain_dict = dict()
    uniprot_chain_dict = dict()
    ligand_chain_dict = dict()
    ligand_name_dict = dict()

    check_m = False
    for line in lines:
        if line[0:10] == 'REMARK 465':
            if line.startswith('REMARK 465   M RES C SSSEQI'):
                check_m = True
                continue
            if check_m:
                res_name3 = line[15:18]
                res_idx = int(line[22:26])
                res_inser = line[22:27]
                chain_id = line[19]
                res_name1 = Res31[res_name3]
                if chain_id not in chain_dict:
                    chain_dict[chain_id] = dict()
                chain_dict[chain_id][res_inser] = res_name1
        elif line[0:6] == 'DBREF ':
            chain_id = line[12]
#            uniprot_code = line[33:41].strip()
            uniprot_id = line[42:52].strip()
            if uniprot_id not in uniprot_chain_dict:
                uniprot_chain_dict[uniprot_id] = list()
            uniprot_chain_dict[uniprot_id] += [chain_id]

        elif line[0:6] == 'MODRES':
            #            res_name3_mod = line[12:15]
            chain_id = line[16]
#            res_idx = int(line[18:22])
            res_inser = line[18:23]
            res_name3 = line[24:27]
            res_name1 = Res31[res_name3]
            if chain_id not in chain_dict:
                chain_dict[chain_id] = dict()
            chain_dict[chain_id][res_inser] = res_name1

        elif line[0:6] == 'HET   ':
            ligand_code = line[7:10]
            chain_id = line[12]
            res_idx = int(line[13:17])
            num_atoms = int(line[20:25])
            if chain_id not in ligand_chain_dict:
                ligand_chain_dict[chain_id] = dict()
            if ligand_code not in ligand_chain_dict[chain_id]:
                ligand_chain_dict[chain_id][ligand_code] = list()
            ligand_chain_dict[chain_id][ligand_code] += [res_idx]

            if ligand_code not in ligand_name_dict:
                ligand_name_dict[ligand_code] = dict()
                ligand_name_dict[ligand_code]['name'] = ''
                ligand_name_dict[ligand_code]['num_atoms'] = num_atoms
        elif line[0:6] == 'HETNAM':
            ligand_code = line[11:14]
            lname = line[15:].strip()
            ligand_name_dict[ligand_code]['name'] += lname

        elif line[0:6] == 'ATOM  ':
            res_name3 = line[17:20]
#            res_idx = int(line[22:26])
            res_inser = line[22:27]

            chain_id = line[21]
            if chain_id not in chain_dict:
                chain_dict[chain_id] = dict()

            res_name1 = Res31[res_name3]
            chain_dict[chain_id][res_inser] = res_name1

    return uniprot_chain_dict, chain_dict, ligand_chain_dict, ligand_name_dict


def read_fasta_uniprot(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    fasta_dict = dict()
    for line in lines:
        if line[0] == '>':
            lis = line.strip().split('|')
            key = lis[1]
            fasta_dict[key] = ''
            uniprot_id = lis[2].strip().split()[0]
        else:
            fasta_dict[key] += line.strip()
    return fasta_dict, uniprot_id


def read_fasta_chain(file_name):
    fp = open(file_name)
    lines = fp.readlines()
    fp.close()

    fasta_dict = dict()
    for line in lines:
        if line[0] == '>':
            lis = line.strip().split('|')
            chains_line = lis[1]

            chains = list()
            idx = chains_line.find('Chains')

            if idx >= 0:
                chain_id = chains_line[idx+7]
                chains += chain_id
                idx = 8
                while True:
                    idx0 = chains_line[idx:].find(',')
                    if idx0 < 0:
                        break

                    chain_id = chains_line[idx+idx0+2]
                    chains += chain_id
                    idx = idx + idx0 + 3
            else:
                idx = chains_line.find('Chain ')
                chain_id = chains_line[idx+6]
                chains += chain_id

            for chain in chains:
                fasta_dict[chain] = ''
        else:
            for chain in chains:
                fasta_dict[chain] += line.strip()
    return fasta_dict


def main():
    if len(sys.argv) < 3:
        print(usage)
        sys.exit()
#    uniprot_code = sys.argv[1]
#    uniprot_code = 'P00533'
#    print(uniprot_code)
    list_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_dir = 'pdb'
    if len(sys.argv) >= 4:
        pdb_dir = sys.argv[3]

#    list_file = '%s_pdb_list.txt' % uniprot_code
#    fasta_file = '%s.fasta' % uniprot_code
    pdb_list = [x.strip().split(';') for x in open(list_file)]

    fasta_dict, uniprot_id = read_fasta_uniprot(fasta_file)
    uniprot_code = list(fasta_dict.keys())[0]
    seq_uniprot = fasta_dict[uniprot_code]

    nw = nwalign.NWalign(g_extend=1.0, g_open2=11.0, lower=1)

    for data in pdb_list:
        pdb_code, method, resolution, chain_ini_fin = data
        line_out = pdb_code
        chain_ids, ini_fin = chain_ini_fin.strip().split('=')

        ini_fin2 = ini_fin.strip().split('-')
        ini = int(ini_fin2[0])
        fin = int(ini_fin2[1])

#        fasta_file = 'pdb/%s.fasta' % (pdb_code)
#        fasta_dict = read_fasta(fasta_file)
        pdb_file = '%s/%s.pdb' % (pdb_dir, pdb_code)
        if not os.path.exists(pdb_file):
            print(pdb_file, 'does not exist')
            continue

        uniprot_chain_dict, chain_dict, ligand_chain_dict, ligand_name_dict = read_pdb(
            pdb_file)

        chain_list = uniprot_chain_dict[uniprot_id]
        chain_out = ''
        for chain_id in chain_list:
            chain_out += '/%s' % chain_id

        line_out += ';%s' % chain_out.strip('/')

        ligand_list = list()
        for chain_id in ligand_chain_dict:
            if chain_id not in chain_list:
                continue
            chain = ligand_chain_dict[chain_id]
            keys = list(chain.keys())
            ligand_list += keys
        ligand_list = list(set(ligand_list))

        ligand_line = ''
        for ligand_id in ligand_list:
            ligand_line += '%s,' % ligand_id
        line_out += ';%s' % ligand_line.strip(',')

        chain_id = chain_list[0]
        res_dict = chain_dict[chain_id]

        res_inser_list = res_dict.keys()
        res_idx_list = list()
        seq_pdb = str()
        check_ins = False
        for resi in res_inser_list:
            res = resi[0:4]
            ins = resi[4]
            if ins != ' ':
                check_ins = True
                break
            res_idx_list += [int(res)]
        if check_ins:
            continue
        res_idx_list = sorted(res_idx_list)
        for idx in res_idx_list:
            resi = '%4d ' % (idx)
            seq_pdb += res_dict[resi]

#        start = min(res_idx_list)
#        end = max(res_idx_list)

        seq2 = seq_uniprot[ini-1:fin]
        gg = len(seq_pdb)-len(seq2)
        seq1 = seq_pdb[gg:]

        check = seq1 == seq2
        if check:
            line_out += ';wild_type'
            print(line_out)
            continue
        align, score = nw.NW(seq_pdb, seq_uniprot)

        q3, d3, q4, d4, count = align[0]
#        print(q3, d3)
        s_count = 0
        check_start = False
        check_end = False
        line_mutation = ''
        gap_count = 0
        for i in range(len(d4)):
            if q4[i] == d4[i]:
                if not check_start:
                    ini = d3[i]
                    check_start = True
                    s_count = 0
                s_count += 1

                check_end = False
                gap_count = 0
                fin = d3[i]
            if q4[i] == '-' and d4[i] != '-':
                if not check_end:
                    check_end = True
                gap_count += 1
            if gap_count >= 10 and s_count <= 5:
                check_start = False
        if not check_end:
            fin = len(seq_uniprot)
        for i in range(ini-1, fin):
            if q4[i] != d4[i]:
                if q4[i].upper() == 'X':
                    continue
                ii = i
                while d3[ii] == '-':
                    ii -= 1
                num = d3[ii]
                line_mutation += '%s%s%s,' % (d4[i], num, q4[i].upper())

        line_mutation = line_mutation.strip(',')
        if line_mutation == '':
            line_mutation = 'wild_type'
        line_out += ';%s' % (line_mutation)

        print(line_out)


if __name__ == "__main__":
    main()
