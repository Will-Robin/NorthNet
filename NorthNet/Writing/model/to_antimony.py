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


def write_inputs_outputs_text(network, input_tokens, output_tokens):
    """
    Write the text for reaction inputs and outputs.

    Parameters
    ----------
    network: NorthNet.Classes.Network
    input_tokens: dict()
    output_tokens: dict()

    Returns
    -------
    io_text: string
    """

    inputs = [input_tokens[i] for i in network.NetworkInputs]
    outputs = [output_tokens[i] for i in network.NetworkOutputs]

    io_text = "\n"
    io_text += "# List of inputs\n"

    io_text += "inputs "
    io_text += " ".join(inputs)
    io_text += ";\n"

    io_text += "\n"
    io_text += "# List of outputs\n"

    io_text += "outputs "
    io_text += " ".join(outputs)
    io_text += ";\n"

    return io_text


def write_reactions_text(
    network,
    compound_tokens,
    rxn_tokens,
    input_tokens,
    input_process_tokens,
    output_tokens,
    output_process_tokens,
    reaction_arrow="=>",
):
    """
    Write the tellurium model reactions lines.

    Parameters
    ----------
    network: list[str]

    compound_tokens: dict

    rxn_tokens: dict

    input_tokens: dict

    output_tokens: dict

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

    reactions_text += "# Reaction inputs\n"
    for c, input in enumerate(network.InputProcesses, c + 1):

        process = network.InputProcesses[input]

        input_source = input_tokens[process.InputID]
        input_compound = [compound_tokens[c] for c in process.InputCompound]

        lhs = input_source
        rhs = ".".join(input_compound)

        equation = f"k{c}*{input_source}"

        input_string = f"{input_process_tokens[input]}: "
        input_string += f"{lhs} {reaction_arrow} {rhs}; "
        input_string += f"{equation}"

        reactions_text += f"{input_string}\n"

    reactions_text += "# Reaction outputs\n"
    for c, output in enumerate(network.OutputProcesses, c + 1):

        process = network.OutputProcesses[output]

        output_source = output_tokens[process.OutputID]
        output_compound = compound_tokens[process.OutputCompound]

        lhs = output_source
        rhs = output_compound

        equation = f"k{c}*{output_compound}"

        output_string = f"{output_process_tokens[output]}: "
        output_string += f"{lhs} {reaction_arrow} {rhs}; "
        output_string += f"{equation}"

        reactions_text += f"{output_string}\n"

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


def write_flow_profile_text(model, input_tokens):
    """
    Write the flow profile for a model in tellurium format.

    Parameters
    ----------
    model: NorthNet.Classes.ModelWriter

    Returns
    -------
    flow_profile_text: str
    """

    flow_profile_text = "\n"

    flow_profile_text += "time: "
    for x in range(0, len(model.flow_profile_time)):
        flow_profile_text += f"{model.flow_profile_time[x]},"

    flow_profile_text = flow_profile_text.strip(",")
    flow_profile_text += "\n\n"

    t0 = model.flow_profile_time[0]
    for x in range(0, len(model.flow_profile_time)):

        t1 = model.flow_profile_time[x]
        flow_profile_text += f"at ({t0} < time < {t1}): "

        t0 = model.flow_profile_time[x]

        for flow in model.flow_profiles:
            flow_token = flow + "_#0"
            input = input_tokens[flow_token]
            flow_value = model.reactor_volume / model.flow_profiles[flow][x]
            flow_profile_text += f"{input}={flow_value},"

        flow_profile_text = flow_profile_text.strip(",")
        flow_profile_text += ";\n"

    flow_profile_text += f"at ({t0} < time): "

    for flow in model.flow_profiles:
        flow_token = flow + "_#0"
        input = input_tokens[flow_token]
        flow_value = model.flow_profiles[flow][x] / model.sigma_flow[x]
        flow_profile_text += f"{input}={flow_value},"

    flow_profile_text = flow_profile_text.strip(",")
    flow_profile_text += ";\n"

    return flow_profile_text


def to_antimony(model, hash_tokens=False):
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
    antimony_text: str
        antimony model.
    """

    reaction_arrow = "=>"
    network = model.network

    compound_tokens = create_token_table(
        network.NetworkCompounds, hash_tokens=hash_tokens
    )

    rxn_tokens = create_token_table(network.NetworkReactions, hash_tokens=hash_tokens)

    input_tokens = create_token_table(network.NetworkInputs, hash_tokens=hash_tokens)
    input_process_tokens = create_token_table(
        network.InputProcesses, hash_tokens=hash_tokens
    )

    experimental_inputs = create_token_table(
        [*model.flow_profiles], hash_tokens=hash_tokens
    )

    for e in experimental_inputs:
        input_tokens[e + "_#0"] = experimental_inputs[e]

    output_tokens = create_token_table(network.NetworkOutputs, hash_tokens=hash_tokens)
    output_process_tokens = create_token_table(
        network.OutputProcesses, hash_tokens=hash_tokens
    )

    compounds = list(compound_tokens.values())

    # List of species
    antimony_text = write_compounds_text(compounds)

    antimony_text += write_inputs_outputs_text(network, input_tokens, output_tokens)

    # List of reactions and reaction rates
    antimony_text += write_reactions_text(
        network,
        compound_tokens,
        rxn_tokens,
        input_tokens,
        input_process_tokens,
        output_tokens,
        output_process_tokens,
        reaction_arrow=reaction_arrow,
    )

    # Initial concentration state variables
    antimony_text += initial_concentrations_text(compounds)

    # Values of kinetic parameters
    antimony_text += write_rate_constant_text(network.NetworkReactions)

    antimony_text += write_flow_profile_text(model, input_tokens)

    return antimony_text
