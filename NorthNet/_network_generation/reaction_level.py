def run_reaction(reactant_compounds, reaction_template):
    '''
    Performs a chemical reaction.

    Mapping can be added:
        Set atom map numbers of reactants before performing the reaction.
        Atom mappings can be removed:
        rdkit.Chem.rdChemReactions.RemoveMappingNumbersFromReactions((ChemicalReaction)reaction) → None

    Parameters
    ----------
    reactants: tuple
        tuple of rdkit Mol objects which take part in the reaction.
    reaction_template: list
        A list of NetGen Reaction_Template objects.

    Returns
    -------
    reactions: list
        A list of NetGen Generated_Reaction objects.
    '''
    reactions = []
    reactants = tuple(r.Mol for r in reactant_compounds)

    reactant_map_numbers = []
    for r in reactants:
        reactant_map_numbers.append([atom.GetAtomMapNum() for atom in r.GetAtoms()])

    ps = reaction_template.Reaction.RunReactants(reactants) # Run the reaction to give a list of products sets
    for u in ps:
        for p in u:
            p = n_gen.incorrect_chiral_H_solve(p)
            Chem.SanitizeMol(p)

        rxn = AllChem.ChemicalReaction() # Create an empty chemical reaction
        [rxn.AddReactantTemplate(r) for r in reactants]
        [rxn.AddProductTemplate(p) for p in u]

        reactions.append( Classes.Generated_Reaction(rxn, reaction_template) )

    return reactions
