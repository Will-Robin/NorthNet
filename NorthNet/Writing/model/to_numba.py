from .numba.topology_to_string import topology_to_string
from .numba.variables_to_string import variables_to_string
from .numba.flow_profile_to_string import flow_profile_to_string


def to_numba(model, numba_decoration=None):
    """
    Convert a model to a numpy/numba code model of a CSTR.

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

    variables_text = variables_to_string(model)

    # Create strings for numba's decoration arguments.
    nf64 = "numba.float64"
    numba_dec = ""
    numba_dec += f"@numba.jit({nf64}[:]({nf64},{nf64}[:],{nf64}[:]),\n"
    numba_dec += "    locals="

    if flow_profile_text == "":
        # No flow profile variables are in the function.
        numba_dec += f"{{'P': {nf64}[:]}}"
    else:
        # Include the flow profile arrays in the decorator.
        numba_dec += f"{{'P': {nf64}[:],'F': {nf64}[:,:],'I':{nf64}[:]}}"

    numba_dec += ",nopython=True)"

    # Create the header import lines
    lines = ["import numpy as np"]
    if numba_decoration == "jit":
        lines.extend(["import numba", "", "", numba_dec])

    elif numba_decoration == "compile":
        lines.extend(["import numba", "from numba.pycc import CC", ""])
        # TODO: add check that the model name is a valid python variable
        mod_name = model.name
        if model.name == "":
            mod_name = "model_func"
        lines.extend([f"cc = CC('{mod_name}')", ""])
        lines.append(
            '@cc.export("model_func", "float64[:](float64,float64[:],float64[:])")'
        )
    else:
        lines.append("")

    # Write the lines for the model function.
    lines.extend(
        ["def model_function(time, S, k):", "", "    P = np.zeros(len(S))", ""]
    )

    # Include flow profiles
    if flow_profile_text != "":
        lines.extend(
            [flow_profile_text, "    i = np.abs(flow_time - time).argmin()", ""]
        )

    # Include the differential equations.
    for m_text in model_text:
        lines.append(f"    {m_text}")

    # Include return value.
    lines.extend(["", "    return P", ""])

    # Lines for the wrapper function.
    lines.append("def wrapper_function(time, S, k):")
    lines.append("    return model_function(time, S, k)")
    lines.append("")

    # Lines for variables
    lines.extend(variables_text)

    # Compilation commands.
    if numba_decoration == "compile":
        lines.extend(["", 'if __name__ == "__main__":', "    cc.compile()"])

    # Convert lines list to a string
    text = "\n".join(lines)

    return text
