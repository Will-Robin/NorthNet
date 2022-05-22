from NorthNet.Utils import utils


def create_token_table(tokens, hash_tokens=False):
    """
    Create key value pairs of tokens as unmodified of as SHA1 hashes.

    Parameters
    ----------
    tokens: list[str]

    Returns
    -------
    token_table: dict
    """

    if hash_tokens:
        token_conversion = lambda x: utils.sha1_hash(x, num_chars=7)
    else:
        token_conversion = lambda x: x

    token_table = {}
    for tok in tokens:
        token_table[tok] = token_conversion(tok)

    return token_table


def write_compounds_text(compounds):
    """
    Write the tellurium model species line.

    Parameters
    ----------
    species: list[str]

    Returns
    -------
    species_text: str
    """

    species_text = "\n"
    species_text += "# List of compounds\n"

    species_text += "species "
    species_text += " ".join(compounds)
    species_text += ";\n"

    return species_text


def write_reactions_text(network, compound_tokens, rxn_tokens, reaction_arrow="=>"):
    """
    Write the tellurium model reactions lines.

    Parameters
    ----------
    network: list[str]

    compound_tokens: dict

    rxn_tokens: dict

    reaction_arrow: string
        "=>" or "->"

    Returns
    -------
    reactions_text: str
    """

    reactions_text = "\n"
    reactions_text += "# List of reactions and rate equations\n"
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

        reactions_text += f"{reaction_string}\n"

    return reactions_text


def initial_concentrations_text(compounds):
    """
    Write the tellurium model compound initialisations.

    Parameters
    ----------
    compounds: list[str]

    Returns
    -------
    conc_text: str
    """

    conc_text = "\n"
    conc_text += "# Initial concentration state variables\n"
    for comp in compounds:
        conc_text += f"{comp} = 0.0;\n"

    return conc_text


def write_rate_constant_text(reactions):
    """
    Write the tellurium model rate constant text.

    Parameters
    ----------
    reactions: list[str]

    Returns
    -------
    k_text: str
    """

    k_text = "\n"
    k_text += "# Values of kinetic parameters\n"
    for c, _ in enumerate(reactions, 1):
        k_text += f"k{c} = 1.0;\n"

    return k_text

def write_flow_profile(model):
    """
    Write tellurium formatted flow profile.

    Parameters
    ----------
    model: Classes.ModelWriter

    Returns
    -------
    flow_text: string
    """
    return ""



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

    compound_tokens = create_token_table(
        network.NetworkCompounds, hash_tokens=hash_tokens
    )

    rxn_tokens = create_token_table(network.NetworkReactions, hash_tokens=hash_tokens)

    compounds = list(compound_tokens.values())

    # List of species
    tellurium_text = write_compounds_text(compounds)

    # List of reactions and reaction rates
    tellurium_text += write_reactions_text(
        network, compound_tokens, rxn_tokens, reaction_arrow=reaction_arrow
    )

    # Initial concentration state variables
    tellurium_text += initial_concentrations_text(compounds)

    # Values of kinetic parameters
    tellurium_text += write_rate_constant_text(network.NetworkReactions)

    return tellurium_text
