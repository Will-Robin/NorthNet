def variables_to_string(model):
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
