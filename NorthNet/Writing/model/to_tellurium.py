def to_tellurium(model):
    """
    Write a model formatted for use with tellurium.

    Parameters
    ----------
    model: Classes.ModelWriter

    Returns
    -------
    tellurium_text: str
    """

    reaction_arrow = "=>"
    network = model.network

    comp_token = lambda x: x
    rxn_token = lambda x: x

    compounds = [comp_token(c) for c in network.NetworkCompounds]

    # List of species
    tellurium_text = "\n"
    tellurium_text += "# List of species\n"

    tellurium_text += "species "
    tellurium_text += " ".join(compounds)
    tellurium_text += ";\n"

    # List of reactions and reaction rates
    tellurium_text += "\n"
    tellurium_text += "# List of reactions and rate equations\n"
    for c,rxn in enumerate(network.NetworkReactions, 1):

        reaction = network.NetworkReactions[rxn]

        reactants = [comp_token(r) for r in reaction.Reactants]
        products  = [comp_token(p) for p in reaction.Products]

        lhs = " + ".join(reactants)
        rhs =  " + ".join(products)

        equation = f"k{c}*" + "*".join(reactants)

        reaction_string = f"{rxn_token(rxn)}: "
        reaction_string += f"{lhs} {reaction_arrow} {rhs}; "
        reaction_string += f"{equation}"

        tellurium_text += f"{reaction_string}\n"

    # Initial concentration state variables
    tellurium_text += "\n"
    tellurium_text += "# Initial concentration state variables\n"
    for comp in compounds:
        tellurium_text += f"{comp_token(comp)} = 0.0;\n"

    # Values of kinetic parameters
    tellurium_text += "\n"
    tellurium_text += "# Values of kinetic parameters\n"
    for c,_ in enumerate(network.NetworkReactions, 1):
        tellurium_text += f"k{c} = 1.0;\n"

    return tellurium_text
