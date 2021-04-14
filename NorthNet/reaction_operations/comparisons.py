path_ReactionDecoderTool = r"C:\Users\willi\Documents\Packages\NorthNet\Reaction_Decoder_Tool"

def compare_reactions(r1,r2, RDT_path = path_ReactionDecoderTool):
    '''
    Uses the Reaction Decoder Tool via system. Figure out what this does later.

    Parameters
    ----------
    r1: str
        Reaction SMILES for first reaction. Format: Reactants>>Products
    r2: str
        Reaction SMILES for second reaction. Format: Reactants>>Products
    RDT_path: str
        Path to Reaction Decoder Tool java executable
    Returns
    -------
    Unfinished
    '''
    os.chdir(RDT_path)
    ### Compare Reactions
    # Compare Reactions using SMILES with precomputed AAM mappings
    run = 'java -jar ReactionDecoder.jar -Q SMI -q "{}"  -T SMI -t "{}" -j COMPARE -f BOTH -u'.format(r1,r2)

    os.system(run)

    # Compare Reactions using RXN files
    'java -jar ReactionDecoder.jar -Q SMI -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH'

    return 0

def recalibrate_atom_environment(favoured_reactions, unfavoured_reactions):
    '''
    Input
    favoured_reactions: NorthNet Reaction Object
    unfavoured_reactions: NorthNet Reaction Object

    Output
    not_modifications: dict {str:[str,]}
    '''
    not_modifications = {}
    for fr in favoured_reactions:
        for ur in unfavoured_reactions:
            #if fr.Reactants != ur.Reactants:
            #    return not_modifications

            fav_reactants   = [Chem.MolFromSmiles(x) for r in favoured_reactions for x in fr.Reactants]
            unfav_reactants = [Chem.MolFromSmiles(x) for r in unfavoured_reactions for x in ur.Reactants]

            fav_products  = [Chem.MolFromSmiles(x) for r in favoured_reactions for x in fr.Products]
            unfav_products  = [Chem.MolFromSmiles(x) for r in unfavoured_reactions for x in ur.Products]

            # Get reacting substructures
            spl = ur.ReactionTemplate.ReactionSMARTS.split('>>')
            reacting_substructs = [Chem.MolFromSmarts(x) for x in spl[0].split('.') if len(Chem.MolFromSmarts(x).GetAtoms()) > 1]
            product_substructs  = [Chem.MolFromSmarts(x) for x in spl[1].split('.') if len(Chem.MolFromSmarts(x).GetAtoms()) > 1]

            # Apply substructure mapping numbers to favoured reactants
            for f in fav_reactants:
                for s in reacting_substructs:
                    molecule_ops.map_substructure_into_mol(f,s)

            # Map favoured products according to substructure labels
            for f in fav_products:
                for s in product_substructs:
                    molecule_ops.map_substructure_into_mol(f,s)

            # Map unfavoured products according to substructure labels
            for u in unfav_products:
                for s in product_substructs:
                    molecule_ops.map_substructure_into_mol(u,s)

            # Get possible spectator substructures
            spectators = []
            for r in fav_reactants:
                for s in reacting_substructs:
                    spectators.extend(molecule_ops.get_reactant_spectators(r, s))

            spectators = list(set(spectators))

            spectators = [Chem.MolFromSmarts(x) for x in spectators]

            # Figure out which spectator group is connected to the reacting atom by searching
            # for it in the reactants. The substructure must be connected to a mapped atom,
            # but not contain any mapped atoms (which would be part of the reaction).
            pos_react_cent_spectators = []
            for s in spectators:
                for u in unfav_products:
                    pos_react_cent_spectators.extend(molecule_ops.get_spectators(u, s))

            # Check through the favoured products and remove spectator groups which occur
            # in both favoured and unfavoured products
            react_cent_spectators = molecule_ops.elimitate_spectators(fav_products, pos_react_cent_spectators)
            react_cent_spectators = pos_react_cent_spectators
            # knowing which groups are attached to reacting atoms in unfavoured products,
            # they can now be used to define unfavoured environments in the reactants.
            reac_substr = [x for x in fr.ReactionTemplate.ReactantSubstructures]

            # dictionary in which to store additions to existing rules by substructure.
            # rules will be added as not recursive SMARTS: !$()
            for x in reac_substr:
                if x not in not_modifications:
                    not_modifications[x] = []

            for r in react_cent_spectators:
                for rc,rs in zip(fav_reactants,reac_substr):
                    for x in molecule_ops.make_new_rule_token(r, rc):
                        not_modifications[rs].append(x)

            return not_modifications
