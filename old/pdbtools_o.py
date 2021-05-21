
    @ classmethod
    def pdbqt_to_pdb(cls, input_file, output_file):
        cls.obabel_rewrite(input_file, output_file)
        model_dict = cls.read_pdb_ligand(output_file)
        model_dict_new = dict()
        keys = model_dict.keys()
        for model_id in keys:
            (pdb_info_lines, ligand_dict, conect_dict) = model_dict[model_id]
            result_new = cls.fix_charge(ligand_dict, conect_dict)
            (ligand_dict_new, conect_dict_new) = result_new
            model_dict_new[model_id] = (
                pdb_info_lines, ligand_dict_new, conect_dict_new)
        cls.write_model_pdb(model_dict_new, output_file)
        cls.obabel_rewrite(output_file, output_file, option='-h')

    @ classmethod
    def fix_ligand(cls, input_file, output_file, neutralize=False, pH=None, is_peptide=False, bond_rebuilding=False, is_fix_atom_idx=False):

        from pymol import cmd
        if bond_rebuilding:
            cls.obabel_rewrite(input_file, output_file, option='-xn')
            cmd.load(output_file)
        else:
            cmd.load(input_file)
        cmd.h_add()
        cmd.save(output_file)

        cls.obabel_rewrite(output_file, output_file)
        model_dict = cls.read_pdb_ligand(output_file)
        model_dict_new = dict()
        keys = model_dict.keys()
        for model_id in keys:
            (pdb_info_lines, ligand_dict, conect_dict) = model_dict[model_id]
            result_new = cls.fix_charge(ligand_dict, conect_dict)
            (ligand_dict_new, conect_dict_new) = result_new
            model_dict_new[model_id] = (
                pdb_info_lines, ligand_dict_new, conect_dict_new)
        cls.write_model_pdb(model_dict_new, output_file)
        if neutralize:
            option1 = '--neutralize'
            cls.obabel_rewrite(output_file, output_file, option=option1)
        if pH is not None:
            option2 = '-p %.1f' % (pH)
            cls.obabel_rewrite(output_file, output_file, option=option2)
            if not is_peptide:
                model_dict = cls.read_pdb_ligand(output_file)
                model_dict_new = dict()
                keys = model_dict.keys()
                for model_id in keys:
                    (pdb_info_lines, ligand_dict,
                     conect_dict) = model_dict[model_id]
                    ligand_dict_new = cls.correct_ligand_pep(ligand_dict)
                    model_dict_new[model_id] = (
                        pdb_info_lines, ligand_dict_new, conect_dict)
                cls.write_model_pdb(model_dict_new, output_file)

        if is_fix_atom_idx and not is_peptide:
            model_dict = cls.read_pdb_ligand(output_file)
            model_dict_new = dict()
            keys = model_dict.keys()
            check_change = False
            for model_id in keys:
                (pdb_info_lines, ligand_dict,
                 conect_dict) = model_dict[model_id]
                ligand_dict_new, check_change0 = cls.fix_pdb_atom_idx(
                    ligand_dict)
                model_dict_new[model_id] = (
                    pdb_info_lines, ligand_dict_new, conect_dict)
                if check_change0:
                    check_change = True
            if check_change:
                cls.write_model_pdb(model_dict_new, output_file)
