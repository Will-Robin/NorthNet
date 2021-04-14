def make_new_rule_token(spectator_tuple, reactant):

    matches = reactant.GetSubstructMatches(spectator_tuple[0])
    # dummy mol will be 'edited' to make new substructure if appropriate
    dummy = deepcopy(reactant)

    for m in matches:
        match_fragment = [dummy.GetAtomWithIdx(i) for i in m]
        intrinsic_maps = [a.GetAtomMapNum() for a in match_fragment] # map numbers in the fragment
        connecting_maps = [n.GetAtomMapNum() for a in match_fragment for n in a.GetNeighbors()] # map numbers around the fragment
        # the found substructure must be attached to the reacting centre,
        # but not contain any atoms involved in the reaction.
        if spectator_tuple[1] in connecting_maps and sum(intrinsic_maps) == 0:
            # the unfavoured environment in the reactant has been found.
            # Create a new environemnt based on it.
            # The rule should consist of the enxisting definition plus
            # a definition of the reacting substructure plus the new
            # environment which does not favour the reaction.
            # Will have to figure out a way of re-writing a recursive
            # SMARTS form.

            # Add atom maps to the dummy as a tag for atoms to keep
            # 22 is easy to see...
            [dummy.GetAtomWithIdx(i).SetAtomMapNum(22) for i in m]
            # Delete all unmapped (untagged) atoms
            [a.SetAtomicNum(0) for a in dummy.GetAtoms() if a.GetAtomMapNum() == 0]
            substruct = Chem.DeleteSubstructs(dummy, Chem.MolFromSmarts('[#0]'))

            # Find the reacting center atom so that the string output cam
            # be rooted at that atom.
            for a in substruct.GetAtoms():
                if a.GetAtomMapNum() == spectator_tuple[1]:
                    center_index_s = a.GetIdx()

            # remove mappings (not needed anymore)
            for a in substruct.GetAtoms():
                if a.GetIdx() == center_index_s:
                    pass
                else:
                    a.SetAtomMapNum(0)

            # The new substructure to be integrated as an exclusion
            new_rule_token = Chem.MolToSmiles(substruct,rootedAtAtom=center_index_s,
                                                canonical=False, allHsExplicit = True,
                                                isomericSmiles = True)
            yield new_rule_token
        else:
            pass

def get_reactant_spectators(molecule, reacting_substructure):
    '''
    Input
    molecule; reacting_substructure: RDKit Mol objects

    Output
    spectators: list of str

    '''
    # Get possible spectator substructures from reactants by deleting all possibly
    # reactive atoms.
    spec_strs = []
    spectators = []
    dummy = deepcopy(molecule) # create copy to 'edit'

    matches = molecule.GetSubstructMatches(reacting_substructure)
    for m in matches:
        [dummy.GetAtomWithIdx(i).SetAtomicNum(0) for i in m]
        edited = Chem.DeleteSubstructs(dummy, Chem.MolFromSmarts('[#0]'))
        frags = Chem.GetMolFrags(edited, asMols=True)
        # a bit of a hack to make sure corrent substructures are generated.
        # should be able to use mols directly, but cannot get substruct
        # matches that way.
        fr_smarts = [Chem.MolToSmarts(f) for f in frags if len(f.GetAtoms()) > 1]
        [spec_strs.append(f) for f in fr_smarts]

    spec_strs = list(set(spec_strs))
    spectators = spec_strs

    return spectators

def get_spectators(molecule, spectator):
    '''
    Input
    molecule: RDKit Mol
    spectator: RDKit Mol

    Output
    react_cent_spectators: list of tuples (Mol, int)

    '''
    # Figure out which spectator group is connected to the reacting atom by searching
    # for it in the reactants. The substructure must be connected to a mapped atom,
    # but not contain anmy mapped atoms (which would be part of the reaction).
    react_cent_spectators = []
    matches = molecule.GetSubstructMatches(spectator)

    for m in matches:
        fragment = [molecule.GetAtomWithIdx(i) for i in m]
        if sum([a.GetAtomMapNum() for a in fragment]) > 0:
            # a mapped atom is present, and the substructure is part of
            # the reaction not a spectator
            pass
        elif sum([n.GetAtomMapNum() for a in fragment for n in a.GetNeighbors()]) > 0:
            # The substructure connects to a reacting atom.
            # Assume only one connection to reacting center is possible
            # (based on current framework, not necessarily true for all
            # chemical reactions).
            # Get the map number of the reacting atom tat  the spectator
            # group is connected to.
            map_num = [n for n in [n.GetAtomMapNum() for a in fragment for n in a.GetNeighbors()] if n > 0]
            # append spectator group and the map number for the atom they
            # are attached to.
            react_cent_spectators.append((spectator,map_num[0]))

        else:
            # it is part of the molecule far from the reacting center
            pass

    return react_cent_spectators

def elimitate_spectators(fav_products, pos_react_cent_spectators):
    react_cent_spectators = []
    for r in pos_react_cent_spectators:
        for f in fav_products:
            matches = f.GetSubstructMatches(r[0])
            if len(matches) == 0:
                react_cent_spectators.append(r)
            else:
                for m in matches:
                    fragment = [f.GetAtomWithIdx(i) for i in m]
                    if sum([a.GetAtomMapNum() for a in fragment]) > 0:
                        # a mapped atom is present, and the substructure is part of
                        # the reaction not a spectator
                        pass
                    elif sum([n.GetAtomMapNum() for a in fragment for n in a.GetNeighbors()]) > 0:
                        # The substructure connects to a reacting atom.
                        # Assume only one connection to reacting center is possible
                        # (based on current framework, not necessarily true for all
                        # chemical reactions).
                        pass
                    else:
                        react_cent_spectators.append(r)

    return react_cent_spectators

def get_substructure_replacements(spectator_group,link_atom_map,reactant):

    matches = reactant.GetSubstructMatches(spectator_group)
    # dummy mol will be 'edited' to make new substructure if appropriate`
    dummy = deepcopy(reactant)

    for m in matches:
        match_fragment = [dummy.GetAtomWithIdx(i) for i in m]
        intrinsic_maps = [a.GetAtomMapNum() for a in match_fragment]
        connecting_maps = [n.GetAtomMapNum() for a in match_fragment for n in a.GetNeighbors()]
        if link_atom_map in connecting_maps and sum(intrinsic_maps) == 0:
            # the unfavoured environment in the reactant has been found.
            # Create a new environemnt based on it.
            # The rule should consist of the enxisting definition plus
            # a definition of the reacting substructure plus the new
            # environment which does not favour the reaction.
            # Will have to figure out a way of re-writing a recursive
            # SMARTS form.

            # Add atom maps to the dummy as a tag for atoms to keep
            # 22 is easy to see...
            [dummy.GetAtomWithIdx(i).SetAtomMapNum(22) for i in m]
            # Delete all unmapped (untagged) atoms
            [a.SetAtomicNum(0) for a in dummy.GetAtoms() if a.GetAtomMapNum() == 0]
            substruct = Chem.DeleteSubstructs(dummy, Chem.MolFromSmarts('[#0]'))

            # now have to figure out how to rewrite all this as a string
            # Find the reacting center atom so that the string output cam
            # be rooted at that atom.
            for a in substruct.GetAtoms():
                if a.GetAtomMapNum() == link_atom_map:
                    center_index_s = a.GetIdx()

            # remove mappings (not needed anymore)
            for a in substruct.GetAtoms():
                a.SetAtomMapNum(0)

            # The new substructure to be integrated as an exclusion
            return Chem.MolToSmiles(substruct,rootedAtAtom=center_index_s,
                                    canonical=False, allHsExplicit = True,
                                    isomericSmiles = True)
        else:
            return None
