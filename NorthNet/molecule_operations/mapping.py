def map_substructure_into_mol(fr,s):
    match_ids = fr.GetSubstructMatches(s)
    for m in match_ids:
        for c,i in enumerate(m):
            mpnum = s.GetAtomWithIdx(c).GetAtomMapNum()
            a = fr.GetAtomWithIdx(i)
            a.SetAtomMapNum(mpnum)

def map_atoms_around_substructure(original, modified, benchmark):
    '''
    original: rdkit Mol object
        Parent substructure

    modified: rdkit Mol object
        Expanded substrcuture
    '''
    modified.UpdatePropertyCache(strict=False) # not sure why this is necessary given what I'm working with, but it prevents rdkit errors

    orig_max_map = max([a.GetAtomMapNum() for a in original.GetAtoms()])

    match_ids = modified.GetSubstructMatches(original)

    for m in match_ids:
        # Set mapping numbers of modified substructure as identical to the
        # original substructure
        for c,i in enumerate(m):
            a = modified.GetAtomWithIdx(i)
            a.SetAtomMapNum(original.GetAtomWithIdx(c).GetAtomMapNum())

        # Remap the rest of the atoms in the modified substructure in order
        # of mapping established above, starting from the highest mappiung number
        # in the original substructure.
        benchmark = molecule_ops.extend_atom_mapping(modified,benchmark = benchmark)

def extend_atom_mapping(mol,benchmark = 0):
    while 0 in [a.GetAtomMapNum() for a in mol.GetAtoms()]:
        molecule_ops.propagate_atom_mapping(mol, benchmark)
        benchmark = max([a.GetAtomMapNum() for a in mol.GetAtoms()])
    return benchmark

def propagate_atom_mapping(mol, benchmark):

    mapped_atoms = [a for a in mol.GetAtoms() if a.HasProp('molAtomMapNumber')]
    mapped_atoms = sorted(mapped_atoms, key = lambda x:x.GetAtomMapNum())
    for a in mapped_atoms:
        for n in sorted([x for x in a.GetNeighbors()], key = lambda x:x.GetAtomicNum()):
            if n.HasProp('molAtomMapNumber'):
                pass
            else:
                benchmark += 1
                n.SetAtomMapNum(benchmark)

    return benchmark

def set_atom_mapping_to_index_one(mol):
    min_mapnum = min([a.GetAtomMapNum() for a in mol.GetAtoms()])

    for a in mol.GetAtoms():
        origMapNum = a.GetAtomMapNum()
        mod = origMapNum - min_mapnum + 1
        a.SetAtomMapNum(mod)

def increment_atom_mapping(mol, increment):
    for a in mol.GetAtoms():
        origMapNum = a.GetAtomMapNum()
        mod = origMapNum + increment
        a.SetAtomMapNum(mod)

def recurs_expand_substruct_within_mol(substruct, molecule, recursion  = 1, benchmark = 0):
    counter = 0
    while counter < recursion:
        if len(substruct.GetAtoms()) == 1:
            # There is no space left to expand.
            break
        substruct = molecule_ops.expand_substruct_within_mol(substruct, molecule, benchmark = benchmark)
        counter += 1

    return substruct

def expand_substruct_within_mol(substruct, molecule, benchmark = 0):

    match_ids = molecule.GetSubstructMatches(substruct)
    links = []

    for m in match_ids:
        for c,i in enumerate(m):
            at = molecule.GetAtomWithIdx(i)
            #at.SetAtomMapNum(substruct.GetAtomWithIdx(c).GetAtomMapNum())
            for b in at.GetBonds():
                links.extend([b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx()])

    links = list(set(links))

    substrate = deepcopy(molecule)
    '''
    mapped_atoms = [a.GetAtomMapNum() for a in substrate.GetAtoms() if a.GetAtomMapNum() > 0]
    if any([mapped_atoms.count(m)>1 for m in mapped_atoms]):
        molecule_ops.discern_mapping_fragment_similarity(substrate, substruct)
    '''
    for a in substrate.GetAtoms():
        if a.GetIdx() in links:
            pass
        else:
            substrate.GetAtomWithIdx(a.GetIdx()).SetAtomicNum(0)

    substrate = Chem.DeleteSubstructs(substrate, Chem.MolFromSmarts('[#0]'))

    return substrate

def discern_mapping_fragment_similarity(mol, substructure):

    mapped_atoms = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum() > 0]
    repeats = set([m for m in mapped_atoms if mapped_atoms.count(m) > 1])
    subs_fp = Chem.RDKFingerprint(substructure)

    bond_idxs = []
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() in repeats:
            for b in a.GetBonds():
                bond_idxs.append(b.GetIdx())

    dummy = deepcopy(mol)
    molecule_ops.add_atom_mapping(dummy, begin_idx = 1)
    best_sim = 0
    best_b_id = bond_idxs[0]
    best_frag = None
    for b in bond_idxs:
        mol1_f = Chem.FragmentOnBonds(dummy, (b,))
        frgs = Chem.GetMolFrags(mol1_f, asMols=True)
        fps = [Chem.RDKFingerprint(x) for x in frgs]
        for c,f in enumerate(fps):
            sim = DataStructs.FingerprintSimilarity(f,subs_fp)
            if sim > best_sim:
                best_sim = sim
                best_b_id = b
                best_frags = frgs

    fps  = [Chem.RDKFingerprint(x) for x in best_frags]
    sims = [DataStructs.FingerprintSimilarity(f,subs_fp) for f in fps]

    remove_map = best_frags[sims.index(min(sims))]

    set_0_idxs = []
    for d in dummy.GetAtoms():
        for r in remove_map.GetAtoms():
            if d.GetAtomMapNum() == r.GetAtomMapNum():
                set_0_idxs.append(d.GetIdx())

    for s in set_0_idxs:
        mol.GetAtomWithIdx(s).SetAtomMapNum(0)

    return 0

def index_map(mol):
    '''
    Cross-reference atom mappings and atom indices.
    Parameters
    ----------
    mol: rdkit Mol object
        Mol object with atom mappings.
    Returns
    -------
    index_map: dict
        A dict [mapping:index]
    '''
    index_map = {}
    for atom in mol.GetAtoms() :
        map_num = atom.GetAtomMapNum()
        if map_num:
            index_map[map_num] = atom.GetIdx()
    return index_map

def add_atom_mapping(mol, begin_idx = 1):
    '''
    Input
    mol: rdkit Mol
    Output
    mol: rdkit Mol
    '''
    for c1,atom in enumerate(mol.GetAtoms(),begin_idx) :
        map_num = atom.SetAtomMapNum(c1)
    return mol

def remove_atom_mapping(mol):
    '''
    Input
    mol: rdkit Mol
    Output
    mol: rdkit Mol
    '''
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol

def fill_empty_map_numbers(reactants):

    max_mapnum = 0
    for r in reactants:
        for a in r.GetAtoms():
            mapnum = a.GetAtomMapNum()
            if mapnum > max_mapnum:
                max_mapnum = mapnum

    max_mapnum += 1
    for r in reactants:
        for a in r.GetAtoms():
            if a.GetAtomMapNum() == 0:
                a.SetAtomMapNum(max_mapnum)
                max_mapnum += 1
