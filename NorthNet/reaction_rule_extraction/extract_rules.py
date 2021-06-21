def extract_reacting_substructures(mapped_reaction_smiles, radius = 1):
    '''
    Get the substructures changed by of a mapped reaction. Prototype.

    Parameters
    ----------
    mapped_reaction_smiles: str
        Reaction SMILES string including mappings.
    radius: int
        radius included in expressing atom environment.

    Returns
    -------
    reaction_string: str
        Reaction SMILES containing only the molecules which changed in enviroment
        across the reaction.
    '''

    changed_atoms = {} # a ditionary in which to store the indices (not mapping) of atoms which change in the reaction for each component.

    prods = [Chem.MolFromSmiles(x) for x in mapped_reaction_smiles.split(">>")[1].split(".")] # creating a list of rdkit.Mol reactants.
    reacs = [Chem.MolFromSmiles(x) for x in mapped_reaction_smiles.split(">>")[0].split(".")] # creating a list of rdkit.Mol products.

    for r in reacs:
        if r not in changed_atoms:
            changed_atoms[r] = [] # Create a list on which to append atom indices.

        for at1 in r.GetAtoms(): # iterate over atoms in reactant molecules.
            env1 = rdmolops.FindAtomEnvironmentOfRadiusN(r,radius,at1.GetIdx(), useHs = False) # Get atom environment
            x = rdmolops.PathToSubmol(r,env1) # convert atom enviroment to string for comparison. Use rdkit test for similarity instead?

            neigh_at1 = [neighbor.GetAtomMapNum() for neighbor in at1.GetNeighbors()]

            for p in prods:
                if p not in changed_atoms:
                    changed_atoms[p] = []

                for at2 in p.GetAtoms():
                    if at1.GetAtomMapNum() == at2.GetAtomMapNum():
                        env2 = rdmolops.FindAtomEnvironmentOfRadiusN(p,radius,at2.GetIdx(), useHs = False)
                        y = rdmolops.PathToSubmol(p,env2)

                        neigh_at2 = [neighbor.GetAtomMapNum() for neighbor in at2.GetNeighbors()]

                        # If the nearest neighbours are different and the atom mappings are the same, the atom has been changed by the reaction.
                        fps = [Chem.RDKFingerprint(m) for m in [x,y]]

                        sim = DataStructs.FingerprintSimilarity(fps[0],fps[1])

                        # Based on neighbouring map numbers changing and based on environment change
                        '''
                        if set(neigh_at1) != set(neigh_at2) or sim < 1.0:
                            changed_atoms[r].append(at1.GetIdx())
                            changed_atoms[p].append(at2.GetIdx())
                        '''
                        # based on environment change
                        if sim < 1.0 or at1.GetDegree() != at2.GetDegree():
                            changed_atoms[r].append(at1.GetIdx())
                            changed_atoms[p].append(at2.GetIdx())

                    else:
                        pass

    reacs_extract = [] # List in which to place substructures which react for reactants.
    for r in reacs:
        z = deepcopy(r) # Copy the molecule so it can be edited safely.
        for atom in z.GetAtoms():
            # The following process deletes atoms from z which do not change in the reaction.
            if atom.GetIdx() in changed_atoms[r]: # If the atom is in the list of changed atoms, ignore it.
                pass
            else:
                # Follows suggestion at https://sourceforge.net/p/rdkit/mailman/message/28157259/
                # Atoms which do not change are changed to dummy atoms.
                z.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(0)
        # Dummy atoms are removed.
        mol3 = Chem.DeleteSubstructs(z, Chem.MolFromSmarts('[#0]'))
        reacs_extract.append(mol3) # Append the SMILES/SMARTS string to the extracted substructures.

    # Perform the same operation as above for products.
    prods_extract = [] # List in which to place substructures which react for products.
    for r in prods:
        z = deepcopy(r)
        for atom in z.GetAtoms():
            if atom.GetIdx() in changed_atoms[r]:
                pass
            else:
                z.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(0)
        mol3 = Chem.DeleteSubstructs(z, Chem.MolFromSmarts('[#0]'))
        prods_extract.append(mol3)

    # Reset atom map numbers

    map_counter = 1
    for mol in reacs_extract:
        for atom in mol.GetAtoms():
            init_mapnum = atom.GetAtomMapNum()
            atom.SetAtomMapNum(map_counter)
            for mol2 in prods_extract:
                for atom2 in mol2.GetAtoms():
                    if atom2.GetAtomMapNum() == init_mapnum:
                        atom2.SetAtomMapNum(map_counter)
            map_counter += 1

    reac_smiles = [Chem.MolToSmiles(r) for r in reacs_extract]
    prod_smiles = [Chem.MolToSmiles(p) for p in prods_extract]

    # Convert the resulting lists into a reaction SMILES/SMARTS string.
    reaction_string = reaction_ops.reconstitute_reaction(reac_smiles,prod_smiles)

    return reaction_string

def extract_reacting_substructures_2(mapped_reaction_smiles, radius = 1):
    '''
    Get the substructures changed by of a mapped reaction. Prototype.

    Parameters
    ----------
    mapped_reaction_smiles: str
        Reaction SMILES string including mappings.
    radius: int
        radius included in expressing atom environment.

    Returns
    -------
    reaction_string: str
        Reaction SMILES containing only the molecules which changed in enviroment
        across the reaction.
    '''

    changed_atoms = {} # a ditionary in which to store the indices (not mapping) of atoms which change in the reaction for each component.

    prods = [Chem.MolFromSmiles(x) for x in mapped_reaction_smiles.split(">>")[1].split(".")] # creating a list of rdkit.Mol reactants.
    reacs = [Chem.MolFromSmiles(x) for x in mapped_reaction_smiles.split(">>")[0].split(".")] # creating a list of rdkit.Mol products.

    for r in reacs:
        if r not in changed_atoms:
            changed_atoms[r] = [] # Create a list on which to append atom indices.

        for at1 in r.GetAtoms(): # iterate over atoms in reactant molecules.
            env1 = rdmolops.FindAtomEnvironmentOfRadiusN(r,radius,at1.GetIdx(), useHs = False) # Get atom environment
            x = rdmolops.PathToSubmol(r,env1) # convert atom enviroment to string for comparison. Use rdkit test for similarity instead?

            neigh_at1 = [neighbor.GetAtomMapNum() for neighbor in at1.GetNeighbors()]

            for p in prods:
                if p not in changed_atoms:
                    changed_atoms[p] = []

                for at2 in p.GetAtoms():
                    if at1.GetAtomMapNum() == at2.GetAtomMapNum():
                        env2 = rdmolops.FindAtomEnvironmentOfRadiusN(p,radius,at2.GetIdx(), useHs = False)
                        y = rdmolops.PathToSubmol(p,env2)

                        neigh_at2 = [neighbor.GetAtomMapNum() for neighbor in at2.GetNeighbors()]

                        # If the nearest neighbours are different and the atom mappings are the same, the atom has been changed by the reaction.
                        fps = [Chem.RDKFingerprint(m) for m in [x,y]]

                        sim = DataStructs.FingerprintSimilarity(fps[0],fps[1])

                        # Based on neighbouring map numbers changing and based on environment change
                        '''
                        if set(neigh_at1) != set(neigh_at2) or sim < 1.0:
                            changed_atoms[r].append(at1.GetIdx())
                            changed_atoms[p].append(at2.GetIdx())
                        '''
                        # based on environment change
                        if sim < 1.0 or at1.GetDegree() != at2.GetDegree():
                            changed_atoms[r].append(at1.GetIdx())
                            changed_atoms[p].append(at2.GetIdx())

                    else:
                        pass

    reacs_extract = [] # List in which to place substructures which react for reactants.
    for r in reacs:
        z = deepcopy(r) # Copy the molecule so it can be edited safely.
        for atom in z.GetAtoms():
            # The following process deletes atoms from z which do not change in the reaction.
            if atom.GetIdx() in changed_atoms[r]: # If the atom is in the list of changed atoms, ignore it.
                pass
            else:
                # Follows suggestion at https://sourceforge.net/p/rdkit/mailman/message/28157259/
                # Atoms which do not change are changed to dummy atoms.
                z.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(0)
        # Dummy atoms are removed.
        mol3 = Chem.DeleteSubstructs(z, Chem.MolFromSmarts('[#0]'))
        reacs_extract.append(mol3) # Append the SMILES/SMARTS string to the extracted substructures.

    # Perform the same operation as above for products.
    prods_extract = [] # List in which to place substructures which react for products.
    for r in prods:
        z = deepcopy(r)
        for atom in z.GetAtoms():
            if atom.GetIdx() in changed_atoms[r]:
                pass
            else:
                z.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(0)
        mol3 = Chem.DeleteSubstructs(z, Chem.MolFromSmarts('[#0]'))
        prods_extract.append(mol3)

    reac_smiles = [Chem.MolToSmiles(r) for r in reacs_extract if Chem.MolToSmiles(r) != '']
    prod_smiles = [Chem.MolToSmiles(p) for p in prods_extract if Chem.MolToSmiles(p) != '']

    # Convert the resulting lists into a reaction SMILES/SMARTS string.
    reaction_string = reaction_ops.reconstitute_reaction(reac_smiles,prod_smiles)

    return reaction_string

def expand_reaction_rule_aam(possible_transformation, mapped_reaction = None, radius = 1):
    '''
    possible_transformation: NorthNet Reaction Object
        An allowed transformation
    radius: int
        Environment to account for
    '''
    if mapped_reaction == None:
        mapped_reaction_smiles = reaction_ops.aam_smiles(possible_transformation.ReactionSMILES)
    else:
        mapped_reaction_smiles = mapped_reaction

    rule_str = reaction_ops.extract_reacting_substructures(mapped_reaction_smiles, radius = radius)

    rxn_template = Classes.Reaction_Template(possible_transformation.ReactionTemplate.Name+'_modified', rule_str, rule_str.split('>>')[0].split("."), rule_str.split('>>')[1].split("."))

    return rxn_template

def expand_reaction_rule(possible_transformation, expansion_steps = 1):
    '''
    possible_transformation: NorthNet Reaction Object
        An allowed transformation
    '''

    rule_str = possible_transformation.ReactionTemplate.ReactionSMARTS
    reactants = [Chem.MolFromSmiles(x) for x in  possible_transformation.Reactants]
    products  = [Chem.MolFromSmiles(x) for x in possible_transformation.Products]

    r_substructs = [Chem.MolFromSmarts(x) for x in rule_str.split('>>')[0].split('.')]
    p_substructs = [Chem.MolFromSmarts(x) for x in rule_str.split('>>')[1].split('.')]


    # assume that the order of substructures is the same as the order of compounds
    # which contain them.
    new_reacs = []
    max_r_mapnum = 0
    max_p_mapnum = 0
    for r in r_substructs:
        max_r_mapnum += len(r.GetAtoms())
    for p in p_substructs:
        max_p_mapnum += len(p.GetAtoms())

    for c1,rs in enumerate(r_substructs):
        rc = reactants[c1]
        new_substr = molecule_ops.recurs_expand_substruct_within_mol(rs, rc, recursion  = expansion_steps)
        #molecule_ops.extend_atom_mapping(new_substr,benchmark = max_r_mapnum)
        new_reacs.append(new_substr)
        max_r_mapnum += (len(new_substr.GetAtoms()) - len(rs.GetAtoms()))
    # assume that the order of substructures is the same as the order of compounds
    # which contain them.
    new_prods = []
    for c1,ps in enumerate(p_substructs):
        pc = products[c1]
        new_substr = molecule_ops.recurs_expand_substruct_within_mol(ps, pc, recursion  = expansion_steps)
        #molecule_ops.extend_atom_mapping(new_substr,benchmark = max_p_mapnum)
        new_prods.append(new_substr)
        max_p_mapnum += (len(new_substr.GetAtoms()) - len(ps.GetAtoms()))

    mapped_reaction = reaction_ops.set_mapping_reactants_products(new_reacs, new_prods, r_substructs, p_substructs)

    r_side = '.'.join([Chem.MolToSmarts(x) for x in new_reacs])
    p_side = '.'.join([Chem.MolToSmarts(x) for x in new_prods])

    rule_str = r_side + '>>' + p_side

    # Hacking to convert atom numbers in smarts to symbols.
    # probably a better way in rdkit, but can't find it yet
    rule_str = rule_str.replace('#6','C')
    rule_str = rule_str.replace('#1','H')
    rule_str = rule_str.replace('#8','O')
    rxn_template = Classes.Reaction_Template(possible_transformation.ReactionTemplate.Name+'_modified', rule_str, rule_str.split('>>')[0].split("."), rule_str.split('>>')[1].split("."))

    return rxn_template

def modify_reaction_rule(reaction_template, replacements):
    '''
    Input
    reaction_template: NorthNet Reaction_Template object
    replacements: dict {str:[str,]}

    Output
    NorthNet Reaction_Template object
    '''
    # Hacking text based replacements
    new_rs = []
    scaffold = reaction_template.ReactionSMARTS.split('>>')[0]
    for r in replacements:
        root_token = Chem.MolFromSmarts(r)
        dummy = Chem.MolFromSmarts(r)
        [a.SetAtomMapNum(0) for a in dummy.GetAtoms()]
        homepos = {a.GetAtomMapNum(): a.GetSmarts() for a in root_token.GetAtoms()}
        token_dict = {a.GetAtomMapNum():["$({})".format(b.GetSmarts())] for a,b in zip(root_token.GetAtoms(), dummy.GetAtoms())}

        if len(replacements[r]) == 0:
            pass
        else:
            # Build a new token
            for x in replacements[r]:
                mol = Chem.MolFromSmarts(x)
                for a in mol.GetAtoms():
                    if a.GetAtomMapNum() > 0:
                        mapnum = a.GetAtomMapNum()
                        a.SetAtomMapNum(0)
                        index = a.GetIdx()
                token_dict[mapnum].append('!$({})'.format(Chem.MolToSmiles(mol,rootedAtAtom=index,canonical=False, allHsExplicit = True)))

            for t in token_dict:
                if len(token_dict[t]) < 2:
                    pass
                else:
                    z = ''.join(token_dict[t])
                    scaffold = scaffold.replace(homepos[t],'[{}:{}]'.format(z,t))

    react_side = scaffold
    prod_side = reaction_template.ReactionSMARTS.split('>>')[1]
    result = '>>'.join([react_side,prod_side])
    from NorthNet import Classes
    return Classes.Reaction_Template(reaction_template.Name+'_modified',result, react_side.split('.'), reaction_template.ProductSubstructures)
