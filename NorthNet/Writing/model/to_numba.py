from .numba.topology_to_string import topology_to_string
from .numba.variables_to_string import variables_to_string
from .numba.flow_profile_to_string import flow_profile_to_string


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

    # Create string representations for the function components
    model_text = topology_to_string(model)

    flow_profile_text = flow_profile_to_string(model, indentation="    ")

    # Create strings for numba's decoration arguments.
    nf64 = "numba.float64"
    numba_dec = ""
    numba_dec += f"@numba.jit({nf64}[:]({nf64},{nf64}[:],{nf64}[:]),\n"
    numba_dec += "\tlocals="

    if flow_profile_text == "":
        numba_dec += f"{{'P': {nf64}[:]}}"
    else:
        numba_dec += f"{{'P': {nf64}[:],'F': {nf64}[:,:],'I':{nf64}[:]}}"

    numba_dec += ",nopython=True)"

    # Create the header import lines
    lines = ["import numpy as np"]
    if numba_decoration == "jit":
        lines.append("import numba\n")
        lines.append("")
        lines.append(numba_dec)

    elif numba_decoration == "compile":
        lines.append("import numba")
        lines.append("from numba.pycc import CC\n")
        # TODO: add check that the model name is a valid python variable
        mod_name = model.name
        if model.name == "":
            mod_name = "model_func"
        lines.append(f"cc = CC('{mod_name}')\n")
        lines.append(
            '@cc.export("model_func", "float64[:](float64,float64[:],float64[:])")'
        )
    else:
        lines.append("")

    # Write the lines for the model function.
    lines.append("def model_function(time, S, k):")
    lines.append("")
    lines.append("    P = np.zeros(len(S))")

    # Include flow profiles
    if flow_profile_text != "":
        lines.append("")
        lines.append(flow_profile_text)
        lines.append("    i = np.abs(flow_time - time).argmin()")
        lines.append("")

    # Include the differential equations.
    for m_text in model_text:
        lines.append(f"    {m_text}")

    # Include return value.
    lines.append("")
    lines.append("    return P")
    lines.append("")
    # Lines for the wrapper function.
    lines.append("def wrapper_function(time, S, k):")
    lines.append("    return model_function(time, S, k)")
    lines.append("")

    # Lines for variables
    lines.extend(variables_to_string(model))

    # Compilation commands.
    if numba_decoration == "compile":
        lines.append("")
        lines.append('if __name__ == "__main__":')
        lines.append("    cc.compile()")

    # Convert lines list to a string
    text = "\n".join(lines)

    return text
