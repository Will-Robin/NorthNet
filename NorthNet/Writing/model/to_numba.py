import numpy as np


def write_flow_profile_text(model, indentation=""):
    """
    Write a series of concentration profiles for each input and a total
    flow rate profile.

    Parameters
    ----------
    model: NorthNet.Classes.ModelWriter
    indentation: str

    Returns
    -------
    text: str
    """

    if len(model.flow_profiles) == 0:
        # No flow profile information
        return ""

    input_concentrations = np.zeros(
        (len(model.flow_profiles), len(model.flow_profile_time))
    )

    for c, flow in enumerate(model.inflows):

        compound = flow.split("_")[0]

        if compound in model.inputs:
            conc = model.inputs[compound]
            flow_rate = model.flow_profiles[flow.split("_")[0]]
            conc_profile = conc * flow_rate / model.reactor_volume
            input_concentrations[c] = conc_profile

    total_flows = np.zeros((len(model.outflows), len(model.sigma_flow)))

    # assume that the output flow rates are equally partitioned between the
    # output channels.
    partitioned_flow = model.sigma_flow / len(model.outflows)
    for c, _ in enumerate(model.outflows):
        total_flows[c] = partitioned_flow / model.reactor_volume

    # Write concentration array
    text = f"{indentation}F_in = np.array(\n"
    text += indentation
    array_text = np.array2string(
        input_concentrations,
        formatter={"float_kind": lambda x: "%.9f" % x},
        separator=",",
        threshold=np.inf,
    )

    text += array_text.replace("\n", f"\n{indentation}")
    text += f"{indentation})\n\n"

    # Write flow profile time axis
    text += f"{indentation}flow_time = np.array(\n"
    text += indentation

    array_text = np.array2string(
        model.flow_profile_time,
        formatter={"float_kind": lambda x: "%.9f" % x},
        separator=",",
        threshold=np.inf,
    )

    text += array_text.replace("\n", f"\n{indentation}")
    text += f"{indentation})\n\n"

    # Write total flow rate.
    text += f"{indentation}total_flow = np.array(\n"
    text += indentation
    array_text = np.array2string(
        total_flows,
        formatter={"float_kind": lambda x: "%.9f" % x},
        separator=",",
        threshold=np.inf,
    )

    text += array_text.replace("\n", f"\n{indentation}")
    text += f"{indentation})\n"

    return text


def write_model_equation_text(model):
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
                input_id = network.InputProcesses[i].InputID
                input_conc = model.inflows[input_id]
                ki = f"+{input_conc}"
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


def write_variables_text(model):
    """
    Write model variables as strings stored in a list

    Parameters
    ----------
    model: NorthNet.Classes.ModelWriter

    Returns
    -------
    lines: list[str]
    """

    get_index = lambda x: int(x[x.find("[") + 1 : x.find("]")])
    lines = []

    lines.append("species = {")
    for k in model.species:
        idx = get_index(model.species[k])
        lines.append(f"'{k}':{idx},")
    lines.append("}")

    lines.append("")
    lines.append("reactions = {")
    for k in model.rate_constants:
        idx = get_index(model.rate_constants[k])
        lines.append(f"'{k}':{idx},")
    lines.append("}")

    lines.append("")
    lines.append("inputs = {")
    for k in model.inputs:
        idx = model.inputs[k]
        lines.append(f"'{k}':{idx},")
    lines.append("}")
    lines.append("")

    lines.append("k = np.zeros(max(reactions.values())+1) # rate constants")
    lines.append("")
    lines.append("S = np.zeros(len(species)) # initial concentrations")

    return lines


def to_numba(model, numba_decoration=None):
    """
    Convert a network to numpy/numba code.

    Parameters
    ----------
    model: NorthNet.Classes.ModelWriter
    numba_decoration: str | None
        Which numba decoration to add.
        "jit": Decoration for just in time compilation.
        "compile": Decoration for ahead of time compilation
        None: No decoration.

    Returns
    -------
    text: str
        The module text.
    """

    flow_profile_text = write_flow_profile_text(model, indentation="    ")

    model_text = write_model_equation_text(model)

    nf64 = "numba.float64"
    numba_dec = ""
    numba_dec += f"@numba.jit({nf64}[:]({nf64},{nf64}[:],{nf64}[:]),\n"
    numba_dec += "\tlocals="

    if flow_profile_text == "":
        numba_dec += f"{{'P': {nf64}[:]}}"
    else:
        numba_dec += f"{{'P': {nf64}[:],'F': {nf64}[:,:],'I':{nf64}[:]}}"

    numba_dec += ",nopython=True)"

    lines = ["import numpy as np"]
    if numba_decoration == "jit":
        lines.append("import numba\n")
        lines.append("")
        lines.append(numba_dec)

    if numba_decoration == "compile":
        lines.append("import numba")
        lines.append("from numba.pycc import CC\n")
        mod_name = model.name
        if model.name == "":
            mod_name = "model_func"
        lines.append(f"cc = CC('{mod_name}')\n")
        lines.append(
            '@cc.export("model_func", "float64[:](float64,float64[:],float64[:])")'
        )

    lines.append("def model_function(time, S, k):")
    lines.append("")
    lines.append("    P = np.zeros(len(S))")

    if flow_profile_text != "":
        lines.append("")
        lines.append(flow_profile_text)
        lines.append("    i = np.abs(flow_time - time).argmin()")
        lines.append("")

    for m_text in model_text:
        lines.append(f"    {m_text}")

    lines.append("")
    lines.append("    return P")
    lines.append("")
    lines.append("def wrapper_function(time, S, k):")
    lines.append("    return model_function(time, S, k)")
    lines.append("")

    lines.extend(write_variables_text(model))

    lines.append("")
    if numba_decoration == "compile":
        lines.append('if __name__ == "__main__":')
        lines.append("    cc.compile()")

    text = "\n".join(lines)

    return text
