from NorthNet import Classes


class Network:
    """
    An object which stored compounds and reactions and the connections
    between them.
    """

    def __init__(self, reactions, name, description):
        """
        The Network object is initialised with a list of Reaction objects. If
        the list is empty, then the network is initialised as an empty network.

        Parameters
        ----------
        reactions: list[NorthNet.Classes.Reaction]
            List of reactions to create the network.
        name: string
            A label name for the network.
        description: string
            A description for the network.

        Attributes
        ----------
        Name: str
            A name for the network.
        Description: str
            A description of the network.
        NetworkReactions: dict
            A dictionary containing the NorthNet.Classes.Reaction objects keyed
            by their reaction SMILES: {reactionSMILES: NorthNet.Classes.Reaction}
        NetworkCompounds: dict
            A dictionary containing the NorthNet.Classes.Compound objects keyed
            by their SMILES: {SMILES: NorthNet.Classes.Compound}
        NetworkInputs: dict
            A dictionary containing the NorthNet.Classes.NetworkInput objects keyed
            by their tokens: {token: NorthNet.Classes.NetworkInput}
        NetworkOutputs: dict
            A dictionary containing the NorthNet.Classes.NetworkOutput objects keyed
            by their tokens: {token: NorthNet.Classes.NetworkOutput}
        InputProcesses: dict
            A dictionary containing the NorthNet.Classes.ReactionInput objects keyed
            by their tokens: {token: NorthNet.Classes.ReactionInput}
        OutputProcesses: dict
            A dictionary containing the NorthNet.Classes.ReactionOutput objects keyed
            by their tokens: {token: NorthNet.Classes.ReactionOutput}
        """

        if isinstance(reactions, list):
            check_rxns = [isinstance(r, Classes.Reaction) for r in reactions]
            assert all(
                check_rxns
            ), """class Network:
                reactions arg should be a list of NorthNet Reaction objects."""
        else:
            assert isinstance(
                reactions, list
            ), """class Network:
                reactions arg should be a list of NorthNet Reaction objects."""

        assert isinstance(
            name, str
        ), """class Network:
                name arg should be a string."""
        assert isinstance(
            name, str
        ), """class Network:
                description arg should be a string."""

        self.Name = name
        self.Description = description

        self.NetworkReactions = dict()

        self.NetworkCompounds = dict()

        self.NetworkInputs = dict()

        self.NetworkOutputs = dict()

        self.InputProcesses = dict()

        self.OutputProcesses = dict()

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_compound(self, compound):
        """
        Add a compound to the network.

        Parameters
        ----------
        compound: NorthNet.Classes.Compound

        Returns
        -------
        None
        """

        assert isinstance(
            compound, Classes.Compound
        ), """Network.add_compound():
                compound arg should be a NorthNet Compound object"""

        if compound.SMILES not in self.NetworkCompounds:
            self.NetworkCompounds[compound.SMILES] = compound

    def add_compounds(self, compounds):
        """
        Add list of NorthNet Compound objects to the Network

        Parameters
        ----------
        compounds: list[NorthNet.Classes.Compound]
            compounds to be added

        Returns
        -------
        None
        """
        if isinstance(compounds, list):
            check_cpds = [isinstance(c, Classes.Compound) for c in compounds]
            assert all(
                check_cpds
            ), """Network.add_compounds():
                    compounds arg should be list of NorthNet Compound objects"""
        else:
            assert isinstance(
                compounds, list
            ), """Network.add_compounds():
                    compounds arg should be list of NorthNet Compound objects"""

        for compound in compounds:
            self.add_compound(compound)

    def remove_compound(self, compound):
        """
        Remove a compound and the reactions in which it is involved
        from the network.

        Parameters
        ----------
        compound: NorthNet.Classes.Compound
            compound to be removed

        Returns
        -------
        None
        """

        assert isinstance(compound, Classes.Compound)

        remove_reactions = []
        remove_reactions.extend(self.NetworkCompounds[compound.SMILES].In)
        remove_reactions.extend(self.NetworkCompounds[compound.SMILES].Out)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions([self.NetworkReactions[r] for r in remove_reactions])

        del self.NetworkCompounds[compound.SMILES]

    def remove_compounds(self, compounds):
        """
        Remove list of compounds and the reactions in which they are involved
        from the network.

        Parameters
        ----------
        compounds: list[NorthNet.Classes.Compound]
            compounds to be removed

        Returns
        -------
        None
        """

        if isinstance(compounds, list):
            check_cpds = [isinstance(c, Classes.Compound) for c in compounds]
            assert all(
                check_cpds
            ), """Network.remove_compounds():
                compounds arg should be a list of NorthNet Compound objects"""
        else:
            assert isinstance(
                compounds, list
            ), """Network.remove_compounds():
                compounds arg should be a list of NorthNet Compound objects"""

        remove_reactions = []
        for compound in compounds:
            remove_reactions.extend(self.NetworkCompounds[compound.SMILES].In)
            remove_reactions.extend(self.NetworkCompounds[compound.SMILES].Out)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions([self.NetworkReactions[r] for r in remove_reactions])

        for compound in compounds:
            del self.NetworkCompounds[compound.SMILES]

    def add_reaction(self, reaction):
        """
        Adds a reaction and its associated reactants and products into the
        Network.

        Parameters
        ----------
        reaction: NorthNet.Classes.Reaction
            reaction to be added

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, Classes.Reaction
        ), """Network.add_reaction:
            reaction arg should be a NorthNet Reaction object."""

        if reaction.ReactionSMILES not in self.NetworkReactions:
            reaction_smiles = reaction.ReactionSMILES

            self.NetworkReactions[reaction_smiles] = reaction

            for reactant in reaction.Reactants:
                if reactant not in self.NetworkCompounds:
                    new_compound = Classes.Compound(reactant)
                    self.add_compound(new_compound)
                    reactant = new_compound.SMILES

                if reaction_smiles not in self.NetworkCompounds[reactant].Out:
                    self.NetworkCompounds[reactant].Out.append(reaction_smiles)

            for product in reaction.Products:
                if product not in self.NetworkCompounds:
                    new_compound = Classes.Compound(product)
                    self.add_compound(new_compound)
                    product = new_compound.SMILES

                if reaction_smiles not in self.NetworkCompounds[product].In:
                    self.NetworkCompounds[product].In.append(reaction_smiles)

    def add_reactions(self, reactions):
        """
        Use the standardised strings information in the Reaction objects
        to build them into the Network.

        Parameters
        ----------
        reactions: list[NorthNet.Classes.Reaction]
            reactions to be added to the Network

        Returns
        -------
        None
        """
        if isinstance(reactions, list):
            check_reactions = [isinstance(c, Classes.Reaction) for c in reactions]
            assert all(
                check_reactions
            ), """Network.add_reactions():
                reactions arg should be a list of NorthNet Reaction objects"""
        else:
            assert isinstance(
                reactions, list
            ), """Network.add_reactions():
                reactions arg should be a list of NorthNet Reaction objects"""

        for reaction in reactions:
            self.add_reaction(reaction)

    def remove_reactions(self, remove_reactions):
        """
        Remove a list of reactions from the Network.

        Parameters
        ----------
        remove_reactions: list[NorthNet.Classes.Reaction] or list[str]
            reactions to be removed (either in reaction SMILES format or as
                NorthNet.Classes.Reaction).

        Returns
        -------
        None
        """

        assertion_msg = """Network.remove_reactions():
            remove_reactions arg should be list of NorthNet Reaction objects"""

        if isinstance(remove_reactions, list):
            check_rxns = [isinstance(c, Classes.Reaction) for c in remove_reactions]
            check_for_strings = [isinstance(r, str) for r in remove_reactions]

            if not all(check_rxns) and not all(check_for_strings):
                assert all(check_rxns), assertion_msg

        else:
            assert isinstance(remove_reactions, list), assertion_msg

        for reaction in remove_reactions:
            if isinstance(reaction, str):
                r_key = reaction
            else:
                r_key = reaction.ReactionSMILES

            for reactant in set(self.NetworkReactions[r_key].Reactants):
                self.NetworkCompounds[reactant].Out.remove(r_key)
            for product in set(self.NetworkReactions[r_key].Products):
                self.NetworkCompounds[product].In.remove(r_key)

            del self.NetworkReactions[r_key]

    def add_input_process(self, input_addition):
        """
        Add a NetworkInput to the Network

        Parameters
        ----------
        input_addition: NorthNet.Classes.ReactionInput
            Input to be added.

        Returns
        -------
        None
        """

        assert isinstance(
            input_addition, Classes.ReactionInput
        ), """Network.add_input_process(): input arg must be a NorthNet
            ReactionInput object"""

        if input_addition.token not in self.InputProcesses:
            for compound in input_addition.InputCompound:
                if compound in self.NetworkCompounds:
                    self.NetworkCompounds[compound].In.append(input_addition.token)

                    self.InputProcesses[input_addition.token] = input_addition

                    if input_addition.InputID not in self.NetworkInputs:
                        self.NetworkInputs[
                            input_addition.InputID
                        ] = Classes.NetworkInput(input_addition.InputID)

                    self.NetworkInputs[input_addition.InputID].Out.append(
                        input_addition.token
                    )

    def add_input_processes(self, inputs):
        """
        For adding ReactionInput to the network

        Parameters
        ----------
        inputs: list[NortNet.Classes.NetworkInput]
            Inputs to be added to the network.

        Returns
        -------
        None
        """

        if isinstance(inputs, list):
            check_inputs = [isinstance(i, Classes.NetworkInput) for i in inputs]

            assert all(
                check_inputs
            ), """Network.add_input_processes(): input arg must be a list of NorthNet
                NetworkInput object"""
        else:
            assert isinstance(
                inputs, list
            ), """Network.a_processesd_inputs(): input arg must be a list of NorthNet
                NetworkInput object"""

        for i in inputs:
            self.add_input_process(i)

    def add_output_process(self, output):
        """
        Add a NetworkOutput to the Network

        Parameters
        ----------
        output NorthNet.Classes.ReactionOutput
            Output to be added

        Returns
        -------
        None
        """
        assert isinstance(
            output, Classes.OutputProcess
        ), """Network.add_output_process(): output arg must be a NorthNet
            ReactionOutput object"""

        if output.OutputCompound in self.NetworkCompounds:

            self.OutputProcesses[output.token] = output

            self.NetworkCompounds[output.OutputCompound].Out.append(output.token)

            if output.OutputID not in self.NetworkOutputs:
                self.NetworkOutputs[output.OutputID] = Classes.NetworkOutput(
                    output.OutputID
                )

            self.NetworkOutputs[output.OutputID].In.append(output.token)

        else:
            # The output compound is not in the network, so cannot be an output
            pass

    def add_output_processes(self, outputs):
        """
        For adding NetworkOutput to the network

        Parameters
        ----------
        outputs: list[NortNet.Classes.NetworkOutput]
            Outputs to be added to the network

        Returns
        -------
        None
        """
        if isinstance(outputs, list):
            check_inputs = [isinstance(i, Classes.NetworkOutput) for i in outputs]
            assert all(
                check_inputs
            ), """Network.add_output_processes(): outputs arg must be a list of NorthNet
                NetworkOutput object"""
        else:
            assert isinstance(
                outputs, list
            ), """Network.add_output_processes():
                outputs arg must be a list of NorthNet NetworkOutput object"""

        for out in outputs:
            self.add_output_process(out)

    def get_reaction(self, reaction):
        """
        Convenience class for getting a reaction using a key

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        NorthNet.Classes.Reaction or None
        """

        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        if reaction in self.NetworkReactions:
            return self.NetworkReactions[reaction]

        print("Reaction not found in Network")
        return None

    def get_reactants(self, reaction):
        """
        Conveniently get the reactants of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.Reactants

    def get_products(self, reaction):
        """
        Conveniently get the products of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.Products

    def get_reaction_template(self, reaction):
        """
        Conveniently get the ReactionTemplate of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.ReactionTemplate

    def get_reaction_SMARTS(self, reaction):
        """
        Conveniently get the Reaction SMARTS of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        reaction_entry = self.get_reaction(reaction)

        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        if reaction_entry.ReactionTemplate is not None:
            return reaction_entry.ReactionTemplate.ReactionSMARTS

        return None

    def get_reaction_name(self, reaction):
        """
        Conveniently get the Name of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions

        Returns
        -------
        None
        """
        assert isinstance(
            reaction, str
        ), """Network.get_reaction(): reaction arg must be a valid
            SMILES string."""

        reaction_entry = self.get_reaction(reaction)

        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        if reaction_entry.ReactionTemplate is not None:
            return reaction_entry.ReactionTemplate.Name

        return None

    def convert_to_networkx(self):
        """
        Converts NorthNet network object to networkx object.

        Parameters
        ----------

        Returns
        -------
        G: networkx.DiGraph
            Networkx version of the NorthNet network.
        """
        import networkx as nx

        Graph = nx.DiGraph()

        for node in self.NetworkCompounds:
            Graph.add_node(node)

        for reaction in self.NetworkReactions:
            for reactant in self.NetworkReactions[reaction].Reactants:
                for product in self.NetworkReactions[reaction].Products:
                    Graph.add_edge(reactant, reaction)
                    Graph.add_edge(reaction, product)

        return Graph
