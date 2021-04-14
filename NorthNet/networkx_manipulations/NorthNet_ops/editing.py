def remove_network_symmetry(network):
    """
    Removes left handed sugars from the network based on the rule:
    if the stereocentre furthest from the carbonyl is S, delete the species and
    its associated reactions from the network. Furthest from the carbonyl is
    defined for now as the last carbon in the canonicalised structure.
    Parameters
    ----------
    Returns
    -------
    None
    """
    node_remove = []
    edge_remove = []

    for n in [*network.NetworkCompounds]:

        ch_centres = Chem.FindMolChiralCenters(network.NetworkCompounds[n].Mol)

        if len(ch_centres) > 0 and ch_centres[-1][-1] == "S":
            node_remove.append(n)
            #print(ch_centres)
            for r in network.NetworkCompounds[n].In:
                edge_remove.append(r)

            for r in network.NetworkCompounds[n].Out:
                edge_remove.append(r)

    node_remove = list(set(node_remove))
    edge_remove = list(set(edge_remove))

    for r in edge_remove:
        for re in network.NetworkReactions[r].Reactants:
            while r in network.NetworkCompounds[re].Out:
                network.NetworkCompounds[re].Out.remove(r)
        for p in network.NetworkReactions[r].Products:
            while r in network.NetworkCompounds[p].In:
                network.NetworkCompounds[p].In.remove(r)

    for r in edge_remove:
        del network.NetworkReactions[r]

    for n in node_remove:
        del network.NetworkCompounds[n]

def add_flow_terms(network, inputs):

    add_reactions = []
    for i in inputs:
        add_reactions.append(">>{}".format(i))
    for c in network.NetworkCompounds:
        add_reactions.append("{}>>".format(c))

    temp_net = n_gen.network_from_reaction_list(add_reactions, network_name = "",
                                        description = "",
                                        reaction_mapping = False)
    flows = []

    for r in temp_net.NetworkReactions:
        flows.append(temp_net.NetworkReactions[r])

    network.add_reactions(flows)

def skip_step(network,substructure):

    e_rem_list = []
    node_rem_list = []
    new_reaction_list = []

    for r in network.NetworkCompounds:


        if network.NetworkCompounds[r].Mol.HasSubstructMatch(substructure.Mol):
            node_rem_list.append(r)
            for i in network.NetworkCompounds[r].In:
                e_rem_list.append(i)

                in_rs = network.NetworkReactions[i].Reactants
                outs = []

                for z in network.NetworkReactions[i].Products:
                    if z.SMILES != r:
                        outs.append(z.SMILES)

                for o in network.NetworkCompounds[r].Out:
                    e_rem_list.append(o)

                    out_ps = network.NetworkReactions[o].Products
                    ins = []

                    for z in network.NetworkReactions[o].Reactants:
                        if z.SMILES != r:
                            ins.append(z.SMILES)

                    LHS = ".".join(ins+in_rs)
                    RHS = ".".join(outs+out_ps)
                    re_str = "{}>>{}".format(LHS,RHS)

                    new_reaction_list.append( Generated_Reaction( Chem.ReactionFromSmiles(re_str), "skip_{}".format(substructure) ) )

    network.add_reactions(new_reaction_list)

    for s in set(e_rem_list):
        del network.NetworkReactions[s]
        for n in network.NetworkCompounds:
            for x in network.NetworkCompounds[n].In:
                if x == s:
                    network.NetworkCompounds[n].In.remove(s)
            for x in network.NetworkCompounds[n].In:
                if x == s:
                    network.NetworkCompounds[n].Out.remove(s)

    for n in set(node_rem_list):
        del network.NetworkCompounds[n]

    return 0

def add_reaction_to_network(network, reactants, reaction_template, exceptions):
    '''
    Extend the network with a specific reaction between specific compounds.

    Parameters
    ----------
    network: NorthNet/NetGen network object
        Network to be extrapolated from.
    reagents: NorthNet/NetGen Compound objects
        Reagents to be applied to the network.
    reaction_template: NetGen Reaction_Template object.
        Reaction template to be used on the network.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''

    reactive_substructs = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]
    reactants = n_gen.reactive_species(reactants, reactive_substructs)
    if len(reactants) > 0:

        if n_gen.check_reaction_input(reactants, reactive_substructs):
            deprot = n_gen.run_reaction(reactants, reaction_template)
            deprot = n_gen.remove_invalid_reactions(deprot, exceptions)

            insertion = [Classes.Reaction(r) for r in deprot]
            network.add_reactions(insertion) # extend the reactions list with the new reactions
        else:
            pass
    else:
        pass
