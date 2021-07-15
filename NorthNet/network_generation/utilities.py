def cleanup(reactions):
    from rdkit import Chem

    '''Used to clear up the output of rdkit reaction function, parsing multiple
    reaction outcomes. Not perfect.
    reactions: list of reaction SMILES strings

    reactions_out: list of cleaned reaction SMILES strings
    '''

    reactions = [r.replace('[H+][O-]','O')    for r in reactions]
    reactions = [r.replace('[O-][H+]','O')    for r in reactions]
    reactions = [r.replace("/"   ,"") for r in reactions]
    reactions = [r.replace( r"\\"  ,"") for r in reactions]


    reactions_out = []
    for c in range(0,len(reactions)):
        # Split the reaction up for inspection
        LHS = reactions[c].split(">>")[0].split(".")
        RHS = reactions[c].split(">>")[1].split(".")
        ins = ''
        for l in LHS:
            mol = Chem.MolFromSmiles(l)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        ins = ins.strip('.')
        ins += '>>'

        for r in RHS:
            mol = Chem.MolFromSmiles(r)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        reactions_out.append(ins.strip('.'))

    return reactions_out

def remove_network_symmetry(network):
    """
    Removes left handed sugars from the network based on the rule:
    if the stereocentre furthest from the carbonyl is S, delete the species and
    its associated reactions from the network. Furthest from the carbonyl is
    defined for now as the last carbon in the canonicalised structure.

    Parameters
    ----------
    network: NorthNet Network object

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
