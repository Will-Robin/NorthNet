import numpy as np
from NorthNet import Classes


class ModelWriter:
    """
    A class designed to generate modelling apparatus by combining a Network
    structure and experimental conditions.
    """

    def __init__(
        self,
        network=None,
        experiment=None,
        conditions=None,
        model_name="",
    ):
        """

        Parameters
        ----------
        network: NorthNet.Classes.Network

        experiment: NorthNet.Classes.DataReport

        conditions: NorthNet.Classes.ExperimentConditions

        Attributes
        ----------
        name: str
        network: NorthNet.Classes.Network
        time: numpy.ndarray[np.float64]
        flow_profile_time: numpy.ndarray[np.float64]
        flow_profiles: dict
        sigma_flow: numpy.ndarray[np.float64]
        reactor_volume: float
        observed_compounds: list[str]
        species: dict
        rate_constants: dict
        inputs: dict
        outputs: dict
        inflows: dict
        outflows: dict
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

        Parameters
        ----------

        Returns
        -------
        None
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
        data_report: NorthNet.Classes.DataReport

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
            input = Classes.InputProcess(f"{inp}_#0>>{inp}")
            self.network.add_input_process(input)

        # No information about multiple outputs for now.
        # Connect a single output to all of the compounds in the network
        for comp in self.network.NetworkCompounds:
            output = Classes.OutputProcess(f"{comp}>>#0")
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
        from NorthNet.Writing import write_flow_profile_text

        text = write_flow_profile_text(self, indentation=indentation)

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

        from NorthNet.Writing import write_model_equation_text

        eq_lines = write_model_equation_text(self)

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

        from NorthNet.Writing import write_variables_text

        lines = write_variables_text(self)

        return lines

    def write_to_module_text(self, numba_decoration=None):
        """

        Parameters
        ----------
        numba_decoration: str or None
            Which numba decoration to add.
            "jit": Decoration for just in time compilation.
            "compile": Decoration for ahead of time compilation
            None: No decoration.

        Returns
        -------
        text: str
            The module text.
        """
        from NorthNet.Writing import model_to_numba

        text = model_to_numba(self, numba_decoration=numba_decoration)

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

    def to_antimony_model(self, token_type="SMILES"):
        """
        Write the model in antimony format.

        Parameters
        ----------
        hash_tokens: bool

        Returns
        -------
        antimony_text: str
        """

        from NorthNet.Writing import model_to_antimony

        antimony_text = model_to_antimony(self, tokens=token_type)

        return antimony_text
