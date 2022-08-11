def topology_to_string(model):
    """
    Writes models as equations with variables which refer to indices of
    arrays:

    P: 1D array of len(model.network.NetworkCompounds)
        Stores the product state if the system following calculation
    S: 1D array of len(model.network.NetworkCompounds)
        Stores the initial state of the system following calculation
        (see model.get_network_tokens())
    k: 1D array of len(model.network.NetworkReactions)
        rate constants arranged in standardised order.
        (see model.get_network_tokens())

    Includes output flow terms for all compounds.

    Parameters
    ----------
    model: NorthNet.Classes.ModelWriter

    Returns
    -------
    eq_lines: list
        List of rate equations in text form.
    """

    network = model.network

    compounds = [*network.NetworkCompounds]

    eq_lines = []

    for count, compound in enumerate(compounds):
        line_text = f"P[{count}] = "
        for i in network.NetworkCompounds[compound].In:
            if i in network.InputProcesses:
                conc = model.inputs[compound]
                input_id = network.InputProcesses[i].InputID
                input_rate = model.inflows[input_id]
                # (flow_rate / vol ) * concentration
                ki = f"+({input_rate}*{conc})"
                line_text += ki
            else:
                reactants = network.NetworkReactions[i].Reactants

                ki = f"+{model.rate_constants[i]}*"

                if len(reactants) == 0:
                    specs = ""
                else:
                    specs = "*".join([model.species[x] for x in reactants])

                line_text += f"{ki}{specs}"

        for out in network.NetworkCompounds[compound].Out:
            if out in network.OutputProcesses:
                output_process = network.OutputProcesses[out]
                out_compound = output_process.OutputCompound
                outlet = output_process.OutputID
                out_flow = model.outflows[outlet]
                reactor_conc = model.species[out_compound]
                ki = f"-{reactor_conc}*{out_flow}"
                line_text += ki
            else:
                ki = f"-{model.rate_constants[out]}*"
                reactants = [
                    model.species[x] for x in network.NetworkReactions[out].Reactants
                ]
                specs = "*".join(reactants)
                line_text += f"{ki}{specs}"

        eq_lines.append(line_text)

    return eq_lines
