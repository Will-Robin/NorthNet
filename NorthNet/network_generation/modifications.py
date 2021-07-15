from NorthNet import Classes

def add_flow_terms(network, inputs):
    '''
    Add flow terms into network (inputs and outputs are
    empty species and )
    '''

    add_inputs = []
    for i in inputs:
        r_obj = Classes.ReactionInput("{}_#0>>{}".format(i,i))
        add_inputs.append(r_obj)

    network.add_inputs(add_inputs)

    add_outputs = []
    for c in network.NetworkCompounds:
        r_obj = Classes.ReactionOutput("{}>>Sample".format(c))
        add_outputs.append(r_obj)

    network.add_outputs(add_outputs)

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
