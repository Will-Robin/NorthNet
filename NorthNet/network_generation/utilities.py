from rdkit import Chem


def cleanup(reactions):
    """
    Used to clear up the output of rdkit reaction function, parsing multiple
    reaction outcomes. Not perfect.

    Parameters
    ----------
    reactions: list[str]
        list of reaction SMILES strings

    Returns
    -------
    reactions_out: list[str]
        list of cleaned reaction SMILES strings
    """

    reactions = [r.replace("[H+][O-]", "O") for r in reactions]
    reactions = [r.replace("[O-][H+]", "O") for r in reactions]
    reactions = [r.replace("/", "") for r in reactions]
    reactions = [r.replace(r"\\", "") for r in reactions]

    reactions_out = []
    for _, reaction in enumerate(reactions):
        # Split the reaction up for inspection
        LHS = reaction.split(">>")[0].split(".")
        RHS = reaction.split(">>")[1].split(".")
        ins = ""
        for l_elem in LHS:
            mol = Chem.MolFromSmiles(l_elem)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            ins += smiles + "."

        ins = ins.strip(".")
        ins += ">>"

        for r_elem in RHS:
            mol = Chem.MolFromSmiles(r_elem)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            ins += smiles + "."

        reactions_out.append(ins.strip("."))

    return reactions_out


def remove_network_symmetry(network):
    """
    Removes left handed sugars from the network based on the rule:
    if the stereocentre furthest from the carbonyl is S, delete the species and
    its associated reactions from the network. Furthest from the carbonyl is
    defined for now as the last carbon in the canonicalised structure.

    Parameters
    ----------
    network: NorthNet.Classes.Network

    Returns
    -------
    None
    """

    node_remove = []
    edge_remove = []

    for compound in [*network.NetworkCompounds]:

        ch_centres = Chem.FindMolChiralCenters(network.NetworkCompounds[compound].Mol)

        if len(ch_centres) > 0 and ch_centres[-1][-1] == "S":
            node_remove.append(compound)

            for reaction in network.NetworkCompounds[compound].In:
                edge_remove.append(reaction)

            for reaction in network.NetworkCompounds[compound].Out:
                edge_remove.append(reaction)

    node_remove = list(set(node_remove))
    edge_remove = list(set(edge_remove))

    for removal in edge_remove:
        for reactant in network.NetworkReactions[removal].Reactants:
            while removal in network.NetworkCompounds[reactant].Out:
                network.NetworkCompounds[reactant].Out.remove(removal)
        for product in network.NetworkReactions[removal].Products:
            while removal in network.NetworkCompounds[product].In:
                network.NetworkCompounds[product].In.remove(removal)

    for removal in edge_remove:
        del network.NetworkReactions[removal]

    for removal in node_remove:
        del network.NetworkCompounds[removal]
