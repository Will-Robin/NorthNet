from rdkit import Chem
from NorthNet import Classes

def add_flow_terms(network, inputs):
    '''
    Add flow terms into network (inputs and outputs are
    empty species and )

    Parameters
    ----------
    network: NorthNet Network object
        Network to add flow terms to.
    inputs: list 
        List of input string tokens to add.    

    Returns
    -------
    None
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
    '''
    Finds compounds in the network which have the specified substructure and
    removes them from the network (including the reactions for which they are
    reactants and products). New reactions between the reactant compounds which
    created the removed compounds and those which are products of reactions of
    the removed compounds are created to fill the gaps.

    Parameters
    ----------
    network: NorthNet Network object
        Network to modify

    Returns
    -------
    None
    '''

    reaction_removal_list = []
    compounds_removal_list = []
    new_reaction_list = []

    for c in network.NetworkCompounds:

        if network.NetworkCompounds[c].Mol.HasSubstructMatch(substructure.Mol):
            compounds_removal_list.append(c)
            # Find reactions connected to the compound to remove
            # and their reactants and products

            # Cycle through the incoming reactions, add them into a list for
            # removal.
            for i in network.NetworkCompounds[c].In:
                reaction_removal_list.append(i)

                # Find the reactants for this reaction and any
                # other products produced with the removed compound
                in_rs = network.NetworkReactions[i].Reactants
                outs = []

                for z in network.NetworkReactions[i].Products:
                    if z.SMILES != c:
                        outs.append(z.SMILES)

                # Go through the outgoing reactions add them into a list for
                # removal
                for o in network.NetworkCompounds[c].Out:
                    reaction_removal_list.append(o)
                    
                    # Find the products for this reaction and any other
                    # reactants which are required with the removed compound.
                    out_ps = network.NetworkReactions[o].Products
                    ins = []

                    for z in network.NetworkReactions[o].Reactants:
                        if z.SMILES != c:
                            ins.append(z.SMILES)

                    # Process the tokens for the new reaction
                    new_reactants = ins+in_rs
                    new_products = outs+out_ps

                    new_reactants.sort()
                    new_products.sort()

                    # Create a new reaction between the incoming reaction's
                    # reactants and the outgoing reaction's products, including
                    # the other compounds required as reactants or products of
                    # both these reactions (which were not tagged for removal
                    LHS = ".".join(new_reactants)
                    RHS = ".".join(new_products)
                    re_str = "{}>>{}".format(LHS,RHS)

                    rdkit_reaction = Chem.ReactionFromSmiles(re_str)
                    new_reaction = Classes.Reaction(rdkit_reaction) 
                    new_reaction_list.append(new_reaction)

    network.add_reactions(new_reaction_list)

    # Remove reactions from the network and any compounds which are exclusively
    # produced by these reactions (if not already picked up in the compounds_removal_list.
    for s in set(reaction_removal_list):
        del network.NetworkReactions[s]
        for n in network.NetworkCompounds:
            for x in network.NetworkCompounds[n].In:
                if x == s:
                    network.NetworkCompounds[n].In.remove(s)
            for x in network.NetworkCompounds[n].In:
                if x == s:
                    network.NetworkCompounds[n].Out.remove(s)

    for n in set(compounds_removal_list):
        del network.NetworkCompounds[n]
