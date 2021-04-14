def get_atom_balance(reaction_object):
    '''
    Parameters
    ----------
    Returns
    -------
    '''
    reactants = reaction_object.ReactionSMILES.split(">>")[0].split(".")
    products = reaction_object.ReactionSMILES.split(">>")[1].split(".")
    species = reactants.extend(products)

    atoms_list = []
    for s in species:
        m = Chem.MolFromSmiles(s)
        for atom in m.GetAtoms():
            atoms_list.append(atom.GetAtomicNum())
    atoms_list = list(set(atoms_list))
    atoms_list = [Chem.MolFromSmarts("[#{}]".format(a)) for a in atoms_list]

    r_count = {Chem.MolToSmiles(a):0 for a in atoms_list}
    p_count = {Chem.MolToSmiles(a):0 for a in atoms_list}
    for a in atoms_list:
        for reac in reactants:
            m = Chem.MolFromSmiles(reac)
            pm = m.GetSubstructMatches(a)
            r_count[Chem.MolToSmiles(a)]+=len(pm)

        for prod in products:
            m = Chem.MolFromSmiles(prod)
            pm = m.GetSubstructMatches(a)
            p_count[Chem.MolToSmiles(a)]+=len(pm)

    return r_count, p_count

def parsing_SMILES_test(reaction_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    tag_list = [] # new tag_list (see above)
    for rx in reaction_list:
        reactants = rx.ReactionSMILES.split(">>")[0].split(".")
        products = rx.ReactionSMILES.split(">>")[1].split(".")
        tag = False
        if None in [Chem.MolFromSmiles(r) for r in reactants]: # if rdkit cannot convert the SMILES string, the string is invalid and None is returned by the function
            tag = True
            rx.Info["Rejection_reason"] = "Rejected: Could not parse SMILES"
        elif None in [Chem.MolFromSmiles(p) for p in products]:
            tag = True
            rx.Info["Rejection_reason"] = "Rejected: Could not parse SMILES"
        tag_list.append(tag)

    test_fails = [x for x,y in zip(reaction_list,tag_list) if y == True]
    reaction_list = [x for x,y in zip(reaction_list,tag_list) if y == False]

    return reaction_list, test_fails

def get_database_entry(reaction_IDs, DB_name = "2019-12-11_Curated_Reaction_Database.txt", DB_storage_dir = r"C:\Users\willi\Documents\PROJECTS\NPC_Project\Databases\Reaxys_Reaction_Database\Processed_Database\2019_12_11"):

        '''
        Parameters
        ----------
        reactions_IDs: list
            List of recation ID strings fow which to search in the database.
        DB_storage_dir: str
            path to directory holding database as tab-delimited text.
        Returns
        -------
        None
        '''
        os.chdir(DB_storage_dir)

        inspect_lines = []

        with open(DB_name, 'r') as f:
            for c,line in enumerate(f):
                if c == 0:
                    ln = line.strip("\n")
                    header = ln.split("\t")
                ln = line.split("\t")
                if ln[0] in reaction_IDs:
                    inspect_lines.append(ln)

        for l in inspect_lines:

            r_SMILES = l[header.index('Reaction')]
            reference = l[header.index('References')]
            print("_______________________________________")
            print("Reaction ID >>", l[header.index('Reaction ID')])
            print("Reaction >>", r_SMILES)
            print("References >>", reference)
            print("Yield >>", l[header.index("Yield (numerical)")])
            print("Reagents >>", l[header.index('Reagent')])
            print("pH >>", l[header.index('pH-Value (Reaction Details)')])
            print("Temperature >>", l[header.index('Temperature (Reaction Details) [C]')])
            print()
            print(l)
            print()
            print("_______________________________________")
