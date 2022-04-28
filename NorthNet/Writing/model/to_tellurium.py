from NorthNet.Utils import utils


def to_tellurium(model, hash_tokens=False):
    """
    Write a model formatted for use with tellurium.

    Parameters
    ----------
    model: Classes.ModelWriter
        Model object containing information.

    hash_tokens: bool
        Whether to use the SHA1 hashes (up to 7 chars) of the reaction network
        tokens.

    Returns
    -------
    tellurium_text: str
        tellurium model.
    """

    reaction_arrow = "=>"
    network = model.network

    compound_tokens = {}
    rxn_tokens = {}
    if hash_tokens:
        for smiles in network.NetworkCompounds:
            compound_tokens[smiles] = utils.sha1_hash(smiles, num_chars=7)
        for reaction in network.NetworkReactions:
            rxn_tokens[reaction] = utils.sha1_hash(reaction, num_chars=7)
    else:
        for smiles in network.NetworkCompounds:
            compound_tokens[smiles] = smiles
        for reaction in network.NetworkReactions:
            rxn_tokens[reaction] = reaction

    compounds = [compound_tokens[c] for c in network.NetworkCompounds]

    # List of species
    tellurium_text = "\n"
    tellurium_text += "# List of species\n"

    tellurium_text += "species "
    tellurium_text += " ".join(compounds)
    tellurium_text += ";\n"

    # List of reactions and reaction rates
    tellurium_text += "\n"
    tellurium_text += "# List of reactions and rate equations\n"
    for c, rxn in enumerate(network.NetworkReactions, 1):

        reaction = network.NetworkReactions[rxn]

        reactants = [compound_tokens[r] for r in reaction.Reactants]
        products = [compound_tokens[p] for p in reaction.Products]

        lhs = " + ".join(reactants)
        rhs = " + ".join(products)

        equation = f"k{c}*" + "*".join(reactants)

        reaction_string = f"{rxn_tokens[rxn]}: "
        reaction_string += f"{lhs} {reaction_arrow} {rhs}; "
        reaction_string += f"{equation}"

        tellurium_text += f"{reaction_string}\n"

    # Initial concentration state variables
    tellurium_text += "\n"
    tellurium_text += "# Initial concentration state variables\n"
    for comp in compounds:
        tellurium_text += f"{comp} = 0.0;\n"

    # Values of kinetic parameters
    tellurium_text += "\n"
    tellurium_text += "# Values of kinetic parameters\n"
    for c, _ in enumerate(network.NetworkReactions, 1):
        tellurium_text += f"k{c} = 1.0;\n"

    return tellurium_text
