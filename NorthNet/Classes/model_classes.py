import sys
import numpy as np
from NorthNet import Classes


class ModelWriter:
    def __init__(
        self,
        network=None,
        experiment=None,
        conditions=None,
        model_name="",
    ):
        """

        A class designed to generate modelling apparatus by combining a Network
        structure and experimental conditions.

        Parameters
        ----------
        network: NorthNet Network

        experiment: NorthNet DataReport

        conditions: NorthNet ExperimentConditions

        """

        if network:
            assert isinstance(
                network, Classes.Network
            ), """Classes.ModelWriter:
                network kwarg should be Network object"""
        if experiment:
            assert isinstance(
                experiment, Classes.DataReport
            ), """Classes.ModelWriter:
                experiment kwarg should be DataReport object."""
        if conditions:
            assert isinstance(
                conditions, Classes.ExperimentConditions
            ), """Classes.ModelWriter:
                conditions kwarg should be ExperimentConditions object."""

        # name
        self.name = model_name

        # Network structure
        self.network = network

        # Reaction conditions details
        self.time = np.array([0.0])
        self.flow_profile_time = np.array([0.0])
        self.flow_profiles = dict()
        self.sigma_flow = []
        self.reactor_volume = 1.0

        # Data attributes
        self.observed_compounds = []

        # Model attributes
        self.species = dict()
        self.rate_constants = dict()
        self.inputs = dict()
        self.outputs = dict()
        self.inflows = dict()
        self.outflows = dict()

        if experiment is None:
            pass
        else:
            self.load_conditions(experiment)

        if conditions is None:
            pass
        else:
            self.load_conditions(conditions)

        if network is None:
            pass
        else:
            self.create_network_tokens()

    def create_network_tokens(self):
        """
        Get dictionaries of tokens for the compounds, reactions, inputs,
        inflows, outflows
        """

        # shorthand for the network
        network = self.network

        # Identies for compounds and reactions
        compounds = [*network.NetworkCompounds]
        reactions = [*network.NetworkReactions]

        # Identies for input and output entities
        input_ids = [*network.NetworkInputs]
        output_ids = [*network.NetworkOutputs]

        # Identies for input and output processes
        network_inputs = [*network.InputProcesses]
        network_outputs = [*network.OutputProcesses]

        # Tokens for reaction concentration terms
        species = {s: f"S[{c}]" for c, s in enumerate(compounds) if s != ""}
        # Tokens for the reaction rate constants of the system
        rate_consts = {k: f"k[{c}]" for c, k in enumerate(reactions)}

        # Tokens for flow terms
        # F_in is an (n x m) array of input concentrations over time
        inflow_rates = {k: f"F_in[{c},i]" for c, k in enumerate(input_ids)}
        outflow_rates = {k: f"total_flow[{c},i]" for c, k in enumerate(output_ids)}

        # Put the data into attributes
        self.species = species
        self.rate_constants = rate_consts
        self.inflows = inflow_rates
        self.outflows = outflow_rates

    def load_conditions(self, data_report):
        """
        Load conditions details into ModelWriter attributes to allow compilation of
        experimental conditions into the model.

        For now, this method will assume that values are in base SI units (e.g.
        M, not mM), and the keys to conditions should follow some relatively strict
        patterns.

        "reactor_volume": gives the reactor volume in L
        "{}/ M": gives the concentration of an inlet in M.
        "{}_flow_{}", not containing "time": gives the flow rate of an input,
        in L/ s
        "{}_flow_time_{}": (contains "flow" and "time") gives the time axis of
        the flow profiles.

        Parameters
        ----------
        conditions: Classes.DataReport

        Returns
        -------
        None
        """

        # Get the series values from the conditions
        self.time = data_report.series_values.copy()

        # Extract flow input information
        for condition in data_report.conditions:
            # Get the reactor volume
            if "reactor_volume" in condition:
                self.reactor_volume = data_report.conditions[condition]

            # Get the input concentrations
            elif "/ M" in condition:
                smiles = condition.split("/")[0]
                self.inputs[smiles] = data_report.conditions[condition]

            # Get the flow profile time
            elif "time" in condition and "flow" in condition:
                self.flow_profile_time = np.array(data_report.conditions[condition])

            # Get flow profile
            elif "flow" in condition and not "time" in condition:
                smiles = condition.split("_")[0]
                self.flow_profiles[smiles] = np.array(data_report.conditions[condition])

        # Get the total flow rate of all the inputs
        self.sigma_flow = np.zeros(len(self.flow_profile_time))
        for flow in self.flow_profiles:
            self.sigma_flow += self.flow_profiles[flow]

        # update the network to reflect the inputs and outputs implied by the
        # conditions
        for inp in self.inputs:
            input = Classes.ReactionInput(f"{inp}_#0>>{inp}")
            self.network.add_input_process(input)

        # No information about multiple outputs for now.
        # Connect a single output to all of the compounds in the network
        for comp in self.network.NetworkCompounds:
            output = Classes.ReactionOutput(f"{comp}>>#0")
            self.network.add_output_process(output)

    def write_flow_profile_text(self, indentation=""):
        """
        Write a series of concentration profiles for each input and a total
        flow rate profile.

        Parameters
        ----------
        indentation: str

        Returns
        -------
        text: str
        """

        if len(self.flow_profiles) == 0:
            # No flow profile information
            return ""

        network = self.network

        input_concentrations = np.zeros(
            (len(self.flow_profiles), len(self.flow_profile_time))
        )

        for c, flow in enumerate(self.inflows):

            compound = flow.split("_")[0]

            if compound in self.inputs:
                conc = self.inputs[compound]
                flow_rate = self.flow_profiles[flow.split("_")[0]]
                conc_profile = conc * flow_rate / self.reactor_volume
                input_concentrations[c] = conc_profile

        total_flows = np.zeros((len(self.outflows), len(self.sigma_flow)))

        # assume that the output flow rates are equally partitioned between the
        # output channels.
        partitioned_flow = self.sigma_flow / len(self.outflows)
        for c, out in enumerate(self.outflows):
            total_flows[c] = partitioned_flow / self.reactor_volume

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
            self.flow_profile_time,
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

    def write_model_equation_text(self):
        """
        Writes models as equations with variables which refer to indices of
        arrays:

        P: 1D array of len(self.network.NetworkCompounds)
            Stores the product state if the system following calculation
        S: 1D array of len(self.network.NetworkCompounds)
            Stores the initial state of the system following calculation
            (see self.get_network_tokens())
        k: 1D array of len(self.network.NetworkReactions)
            rate constants arranged in standardised order.
            (see self.get_network_tokens())

        Includes output flow terms for all compounds.

        Parameters
        ----------
        None

        Returns
        -------
        eq_lines: list
            List of rate equations in text form.
        """

        network = self.network

        compounds = [*network.NetworkCompounds]

        eq_lines = []

        for count, compound in enumerate(compounds):
            line_text = f"P[{count}] = "
            for i in network.NetworkCompounds[compound].In:
                if i in network.InputProcesses:
                    input_id = network.InputProcesses[i].InputID
                    input_conc = self.inflows[input_id]
                    ki = f"+{input_conc}"
                    line_text += ki
                else:
                    reactants = network.NetworkReactions[i].Reactants

                    ki = f"+{self.rate_constants[i]}*"

                    if len(reactants) == 0:
                        specs = ""
                    else:
                        specs = "*".join([self.species[x] for x in reactants])

                    line_text += f"{ki}{specs}"

            for out in network.NetworkCompounds[compound].Out:
                if out in network.OutputProcesses:
                    output_process = network.OutputProcesses[out]
                    out_compound = output_process.OutputCompound
                    outlet = output_process.OutputID
                    out_flow = self.outflows[outlet]
                    reactor_conc = self.species[out_compound]
                    ki = f"-{reactor_conc}*{out_flow}"
                    line_text += ki
                else:
                    ki = f"-{self.rate_constants[out]}*"
                    reactants = [
                        self.species[x] for x in network.NetworkReactions[out].Reactants
                    ]
                    specs = "*".join(reactants)
                    line_text += f"{ki}{specs}"

            eq_lines.append(line_text)

        return eq_lines

    def write_variables_text(self):
        """
        Write model variables as strings stored in a list

        Parameters
        ----------
        None

        Returns
        -------
        lines: list[str]

        """

        get_index = lambda x: int(x[x.find("[") + 1 : x.find("]")])
        lines = []

        lines.append("species = {")
        for k in self.species:
            idx = get_index(self.species[k])
            lines.append(f"'{k}':{idx},")
        lines.append("}")

        lines.append("")
        lines.append("reactions = {")
        for k in self.rate_constants:
            idx = get_index(self.rate_constants[k])
            lines.append(f"'{k}':{idx},")
        lines.append("}")

        lines.append("")
        lines.append("inputs = {")
        for k in self.inputs:
            idx = self.inputs[k]
            lines.append(f"'{k}':{idx},")
        lines.append("}")
        lines.append("")

        lines.append("k = np.zeros(max(reactions.values())+1) # rate constants")
        lines.append("")
        lines.append("S = np.zeros(len(species)) # initial concentrations")

        return lines

    def write_to_module_text(self, numba_decoration=None):
        """

        Parameters
        ----------
        numba_decoration: str
            Which numba decoration to add.
            "jit": Decoration for just in time compilation.
            "compile": Decoration for ahead of time compilation
            None: No decoration.

        Returns
        -------
        text: str
            The module text.
        """

        flow_profile_text = self.write_flow_profile_text(indentation="    ")

        model_text = self.write_model_equation_text()

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
            mod_name = self.name
            if self.name == "":
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

        lines.append("    P *= time")
        lines.append("")
        lines.append("    return P")
        lines.append("")
        lines.append("def wrapper_function(time, S, k):")
        lines.append("    return model_function(time, S, k)")
        lines.append("")

        lines.extend(self.write_variables_text())

        lines.append("")
        if numba_decoration == "compile":
            lines.append('if __name__ == "__main__":')
            lines.append("    cc.compile()")

        text = "\n".join(lines)

        return text

    def write_model_matrix_text(self):
        """

        Prototype for writing the model as an array.

        Parameters
        ----------

        Returns
        -------
        mat_text: str
            Rate equations in numpy matrix form.
        """
        compounds = [*self.network.NetworkCompounds]

        species = self.species
        rate_consts = self.rate_constants
        inflows = self.inputs

        ratemat = [["0" for _ in species] for _ in species]

        for compound in compounds:
            ind1 = compounds.index(compound)

            """outgoing reactions"""
            for i in self.network.NetworkCompounds[compound].In:

                token = ""
                ind2 = ind1
                reacs = self.network.NetworkReactions[i].Reactants
                ki = rate_consts[i]

                if len(reacs) == 0:
                    token = f"(+{ki}*{inflows[compound]})/{species[compound]}"
                    ind2 = ind1

                if len(reacs) == 1:
                    n2 = reacs[0]
                    ind2 = compounds.index(n2)
                    token = "+" + ki

                if len(reacs) == 2:
                    n2 = reacs[0]
                    n3 = reacs[1]
                    ind2 = compounds.index(n3)
                    token = "+" + ki + "*" + species[n2]

                if len(reacs) == 3:
                    n2 = reacs[0]
                    n3 = reacs[1]
                    n4 = reacs[2]
                    ind2 = compounds.index(n3)
                    token = "+" + ki + "*" + species[n2] + "*" + species[n4]

                ratemat[ind1][ind2] += token

            for out in self.network.NetworkCompounds[compound].Out:

                token = ""
                ind2 = ind1
                reacs = self.network.NetworkReactions[out].Reactants
                ki = rate_consts[out]

                if len(reacs) == 1:
                    token = "-" + ki
                    ind2 = ind1

                if len(reacs) == 2:
                    z = reacs[:]
                    z.remove(compound)
                    n2 = z[0]
                    ind2 = compounds.index(n2)
                    token = "-" + ki + "*" + species[compound]

                if len(reacs) == 3:
                    z = reacs[:]
                    z.remove(compound)
                    n2 = z[0]
                    n3 = z[1]
                    ind2 = compounds.index(n2)
                    token = "-" + ki + "*" + species[compound] + "*" + species[n3]

                ratemat[ind1][ind2] += token

        mat_text = "["
        for element in ratemat:
            mat_text += "[" + ",".join(element) + "],\n"

        mat_text = mat_text.strip(",\n") + "]"

        return mat_text

    def to_tellurium_model(self, hash_tokens=False):
        """
        Write the model in a format for use with tellurium.
        """

        from NorthNet.Writing import model_to_tellurium

        tellurium_text = model_to_tellurium(self, hash_tokens=hash_tokens)

        return tellurium_text
