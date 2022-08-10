import numpy as np


def flow_profile_to_string(model, indentation=""):
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
