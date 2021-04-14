path_ReactionDecoderTool = r"C:\Users\willi\Documents\Packages\NorthNet\Reaction_Decoder_Tool"

def aam_smiles(reaction_smiles, RDT_path = path_ReactionDecoderTool):
    '''
    Maps a SMILES reaction by calling the Reaction Decoder Tool via system.

    Parameters
    ----------
    reaction_smiles: str
        Reaction SMILES for reaction to be mapped. Format: Reactants>>Products
    RDT_path: str
        Path to Reaction Decoder Tool java executable
    Returns
    -------
    mapped_reaction: str
        Mapped reaction SMILES.
    '''
    wd = os.getcwd()
    os.chdir(RDT_path)

    #AMM using SMILES
    aam_smiles = 'java -jar ReactionDecoder.jar -Q SMI -q "{}" -g -j AAM -f TEXT'.format(reaction_smiles)

    os.system(aam_smiles)

    with open("ECBLAST_smiles_AAM.txt", "r") as f:
        for line in f:
            if "SELECTED AAM MAPPING" in line:
                line = next(f)
                mapped_reaction = line.strip("\n")
                break

    os.chdir(wd)

    return mapped_reaction

def annotate_smiles(reaction_smiles, RDT_path = path_ReactionDecoderTool):
    '''
    Annotate a reaction given in SMILES format by calling the Reaction Decoder Tool via system.

    Parameters
    ----------
    reaction_smiles: str
        Reaction SMILES for reaction to be annotated. Format: Reactants>>Products
    RDT_path: str
        Path to Reaction Decoder Tool java executable
    Returns
    -------
    Unfinished.
    '''
    os.chdir(RDT_path)
    ### Annotate Reaction using SMILES
    ann_reac_smiles = 'java -jar ReactionDecoder.jar -Q SMI -q "{}" -g -j ANNOTATE -f XML'.format(reaction_smiles)

    os.system(ann_reac_smiles)

    return 0

def set_mapping_reactants_products(reactants, products, reactant_substructures, product_substructures):

    # Match substructures in products and reactants
    '''
    Finding non-reacting groups attached to core reacting functional groups
    using fragmentation.
    '''
    r_bond_breaks = []
    max_mapnum = 0
    for r,s in zip(reactants, reactant_substructures):
        match_ids = r.GetSubstructMatches(s)
        bidx_list = []
        for m in match_ids:
            for c,i in enumerate(m):
                at = r.GetAtomWithIdx(i)
                mapnum = s.GetAtomWithIdx(c).GetAtomMapNum()
                at.SetAtomMapNum(mapnum)
                bidx_list.extend([b.GetIdx() for b in at.GetBonds()])


        bidx_list = list(set(bidx_list))
        r_bond_breaks.append(bidx_list)

    supp_react_frags = []
    for c,bonds in enumerate(r_bond_breaks):
        r = reactants[c]
        s = reactant_substructures[c]

        subs_fp = Chem.RDKFingerprint(s)

        best_sim = 0
        best_b_id = 0
        best_frag = None

        for b in bonds:
            mol1_f = Chem.FragmentOnBonds(r, (b,))
            frgs = Chem.GetMolFrags(mol1_f, asMols=True)
            fps = [Chem.RDKFingerprint(x) for x in frgs]
            for c,f in enumerate(fps):
                sim = DataStructs.FingerprintSimilarity(f,subs_fp)
                if sim > best_sim:
                    best_sim = sim
                    best_b_id = b
                    best_frags = frgs

        [supp_react_frags.append(x) for x in best_frags if not x.HasSubstructMatch(s)]

    # Give atom map numbers to unmapped reactant atoms
    #molecule_ops.fill_empty_map_numbers(reactants)

    # Create a log of neighbouring atoms by map number.
    map_neighbours = {}
    for r in reactants:
        for a in r.GetAtoms():
            mapnum = a.GetAtomMapNum()
            map_neighbours[mapnum] = [n.GetAtomMapNum() for n in a.GetNeighbors()]

    '''
    Finding non-reacting groups in the products via comparison with fragments
    with the expected reacting substructure and the non-reacting groups identified
    from the reactants.
    '''
    p_bond_breaks = []
    for p,s in zip(products, product_substructures):
        match_ids = p.GetSubstructMatches(s)
        bidx_list = []
        for m in match_ids:
            for c,i in enumerate(m):
                at = p.GetAtomWithIdx(i)
                bidx_list.extend([b.GetIdx() for b in at.GetBonds()])

        bidx_list = list(set(bidx_list))
        p_bond_breaks.append(bidx_list)

    # create dummy products with its sown atom mapping so that atoms can be
    # followed through fragmentations.
    dummy_products = [deepcopy(p) for p in products]
    supp_prod_frags = []
    for c,bonds in enumerate(p_bond_breaks):
        p = dummy_products[c]
        molecule_ops.add_atom_mapping(p, begin_idx = 1)

        subs_fp = Chem.RDKFingerprint(s)

        best_sim = 0
        best_b_id = 0
        best_frag = None

        test_frags = supp_react_frags + [product_substructures[c]]

        for b1 in bonds:
            for b2 in bonds:
                if b1 == b2:
                    continue
                mol1_f = Chem.FragmentOnBonds(p, (b1,b2))
                frgs = Chem.GetMolFrags(mol1_f, asMols=True)
                fps = [Chem.RDKFingerprint(x) for x in frgs]
                fp_score = 0

                for rs in test_frags:
                    subs_fp = Chem.RDKFingerprint(rs)
                    for f in fps:
                        sim = DataStructs.FingerprintSimilarity(f,subs_fp)
                        fp_score += sim
                    if fp_score > best_sim:
                        best_sim = fp_score
                        best_b_id = (b1,b2)
                        best_frags = frgs

        [supp_prod_frags.append(x) for x in best_frags]

    '''
    Apply mappings based on the identified fragments from the previous two steps
    '''
    # First find reactive substrcuture
    for c1, s in enumerate(product_substructures):
        for p in supp_prod_frags:
            match_ids = p.GetSubstructMatches(s)
            for m in match_ids:
                for c2,i in enumerate(m):
                    a = s.GetAtomWithIdx(c2)
                    transfer_mapnum = a.GetAtomMapNum()

                    ref_mapnum = p.GetAtomWithIdx(i).GetAtomMapNum()
                    for refat in dummy_products[c1].GetAtoms():
                        if refat.GetAtomMapNum() == ref_mapnum:
                            refidx = refat.GetIdx()
                            products[c1].GetAtomWithIdx(refidx).SetAtomMapNum(transfer_mapnum)
                            break
    '''
    # Assuming bonds to unmapped atoms are left intact, map according to
    # relative arrangements in reactants
    for p in products:

        mapped = [a.GetAtomMapNum() for a in p.GetAtoms()]

        for a in p.GetAtoms():
            if a.GetAtomMapNum() == 0:
                neigh_maps = [n.GetAtomMapNum() for n in a.GetNeighbors() if n.GetAtomMapNum() > 0]
                possible_neighbours = [item for n in neigh_maps for item in map_neighbours[n]]
                for pn in possible_neighbours:
                    if pn in mapped:
                        pass
                    else:
                        a.SetAtomMapNum(pn)
    '''
    return 0
