import networkx as nx
from NorthNet import Classes


class SubstructureNetwork:
    """
    A network designed to show the relationship between functional group
    transformations.

    The network has three kinds of node: reaction rules, compounds, and
    substructures. Substructures connect reactions to compounds. Compounds and
    reactions can only link to substructures.
    """

    def __init__(self, reactions, name, description):
        """
        Parameters
        ----------
        reactions: list[NorthNet.Classes.Reaction]
            List of reactions containing NorthNet.Classes.ReactionTemplates.
        name: str
            A name for the network.
        description: str
            A description of the network.

        Attributes
        ----------
        Name: str
            A name for the network.
        Description: str
            A description of the network.
        Substructures: dict
            Dictionary containing NorthNet.Classes.Substructure objects keyed
            by SMARTS.
        ReactionRules: dict
            Dictionary containing NorthNet.Classes.ReactionTemplate objects keyed
            by reaction SMARTS.
        Compounds: dict
            Dictionary containing NorthNet.Classes.Compound objects keyed
            by SMILES.
        """

        assert isinstance(
            name, str
        ), """class SubstructureNetwork:
            name arg should be a string."""
        assert isinstance(
            name, str
        ), """class SubstructureNetwork:
            description arg should be a string."""

        self.Name = name
        self.Description = description
        self.Substructures = dict()
        self.ReactionRules = dict()
        self.Compounds = dict()

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_compound(self, compound):
        """
        Add a compound to the SubstructureNetwork.

        The compound will only be added to the SubstructureNetwork if it has at
        least one substructure match to one of the substructures already in the
        network.

        Parameters
        ----------
        compound: NorthNet.Classes.Compound

        Returns
        -------
        None
        """

        assert isinstance(
            compound, Classes.Compound
        ), """SubstructureNetwork.add_compound():
                compound arg should be a NorthNet Compound object"""

        # Only perform insertion if the compound is not already in the
        # SubstructureNetwork
        if compound.SMILES not in self.Compounds:
            # find the subtructures in the SubstructureNetwork which match the
            # compound

            add_to_net = False  # changes to True if a match is found
            for substruct in self.Substructures:
                substructure = self.Substructures[substruct].Mol
                if compound.Mol.HasSubstructMatch(substructure):
                    compound.ReactiveSubstructures.append(substruct)
                    substructure.MatchingCompounds.append(compound.SMILES)
                    if not add_to_net:
                        add_to_net = True

            if add_to_net:
                self.Compounds[compound.SMILES] = compound

    def add_compounds(self, compounds):
        """
        Add list of NorthNet Compound objects to the SubstructureNetwork

        A compound will only be added to the SubstructureNetwork if it has at
        least one substructure match to one of the substructures already in the
        network (see self.add_compound()).

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
            ), """SNetwork.add_compounds():
                    compounds arg should be list of NorthNet Compound objects"""
        else:
            assert isinstance(
                compounds, list
            ), """SNetwork.add_compounds():
                    compounds arg should be list of NorthNet Compound objects"""

        for compound in compounds:
            self.add_compound(compound)

    def remove_compound(self, compound):
        """
        Remove a compound from the substructure network.


        Parameters
        ----------
        compound: NorthNet.Classes.Compound

        Returns
        -------
        None
        """

        assert isinstance(
            compound, Classes.Compound
        ), """SubstructureNetwork.remove_compound:
            compound arg should be a NorthNet Compound object."""

        # remove the compound connection to substructures
        for substruct in compound.ReactiveSubstructures:
            self.Substructures[substruct].MatchingCompounds.remove(compound)

        del self.Compounds[compound.SMILES]

    def remove_compounds(self, compounds):
        """
        Remove list of compounds from the SubstructureNetwork.

        Parameters
        ----------
        compounds: list[NorthNet.Classes.Compound]
            compounds to be removed
        garbage_collection: bool
            Whether to run 'garbage collection' to remove unconnected nodes
            from the network.

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

        for compound in compounds:
            self.remove_compound(compound)

    def add_reaction_rule(self, rule):
        """
        Add a reaction rule to the network.

        Parameters
        ----------
        rule: NorthNet.Classes.ReactionTemplate

        Returns
        -------
        None
        """

        assert isinstance(
            rule, Classes.ReactionTemplate
        ), """SubstructureNetwork.add_template:
            rule arg should be a NorthNet ReactionTemplate object."""

        rxn_key = rule.ReactionSMARTS

        for rxn_substr in rule.ReactantSubstructures:

            # connect to reaction
            if rxn_substr in self.Substructures:
                pass
            else:
                self.Substructures[rxn_substr] = Classes.Substructure(rxn_substr)

            self.Substructures[rxn_substr].ReactionParticipations.append(rxn_key)

            working_substruct = self.Substructures[rxn_substr].Mol

            # connect to compounds
            for comp in self.Compounds:
                compound = self.Compounds[comp].Mol
                if compound.HasSubstructMatch(working_substruct):
                    self.Compounds[comp].ReactiveSubstructures.append(rxn_substr)
                    self.Substructures[rxn_substr].MatchingCompounds.append(comp)

        for p_substruct in rule.ProductSubstructures:

            # connect to reaction
            if p_substruct in self.Substructures:
                pass
            else:
                self.Substructures[p_substruct] = Classes.Substructure(p_substruct)

            self.Substructures[p_substruct].ProducingReactions.append(rxn_key)

            working_substruct = self.Substructures[p_substruct].Mol

            # connect to compounds
            for comp in self.Compounds:
                compound = self.Compounds[comp].Mol
                if compound.HasSubstructMatch(working_substruct):
                    self.Compounds[comp].ReactiveSubstructures.append(p_substruct)
                    self.Substructures[p_substruct].MatchingCompounds.append(comp)

    def add_reaction_rules(self, rules):
        """
        Add reaction templates to the network.

        Parameters
        ----------
        rules: list[NorthNet.Classes.ReactionTemplate]

        Returns
        -------
        None
        """

        if isinstance(rules, list):
            check_templates = [isinstance(c, Classes.ReactionTemplate) for c in rules]
            assert all(
                check_templates
            ), """SubstructureNetwork.add_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""
        else:
            assert isinstance(
                rules, list
            ), """SubstructureNetwork.add_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""

        for t in rules:
            if t.ReactionSMARTS in self.ReactionRules:
                pass
            else:
                self.add_reaction_rule(t)

    def remove_reaction_rule(self, rule):
        """
        Remove a reaction rule from a SubstructureNetwork.

        Parameters
        ----------
        rule: NorthNet.Classes.ReactionTemplate

        Returns
        -------
        None
        """

        assert isinstance(
            rule, Classes.ReactionTemplate
        ), """SubstructureNetwork.remove_template:
            template arg should be a NorthNet ReactionTemplate object."""

        # Disconnect the reaction rule from substructures
        for substr in rule.ReactantSubstructures:
            substructure = self.Substructures[substr]
            rxn_partic = substructure.ReactionParticipations
            rxn_partic.remove(rule.ReactionSMARTS)

        for substr in rule.ReactantSubstructures:
            substructure = self.Substructures[substr]
            prod_rxns = substructure.ProducingReactions
            prod_rxns.remove(rule.ReactionSMARTS)

        # Check for disconnected substructures
        remove_substructs = []
        for substr in self.Substructures:
            substructure = self.Substructures[substr]
            rxn_partic = substructure.ReactionParticipations
            prod_rxns = substructure.ProducingReactions
            if len(prod_rxns) == 0 and len(rxn_partic) == 0:
                remove_substructs.append(substr)

        # disconnect the substructure from any matching compounds.
        for r_substr in remove_substructs:
            for comp in self.Substructures[r_substr].MatchingCompounds:
                self.Compounds[comp].ReactiveSubstructures.remove(r_substr)
            del self.Substructures[r_substr]

        # Check for disconnected compounds
        remove_compounds = []
        for comp in self.Compounds:
            compound = self.Compounds[comp]
            if len(compound.ReactiveSubstructures) == 0:
                remove_compounds.append(comp)

        for comp in remove_compounds:
            del self.Compounds[comp]

        # Finally, remove the reaction rule
        del self.ReactionRules[rule.ReactionSMARTS]

    def remove_reaction_rules(self, rules):
        """
        Remove a list of ReactionTemplates from the SubstructureNetwork.

        Parameters
        ----------
        rules: list[NorthNet.Classes.ReactionTemplate]

        Returns
        -------
        None
        """

        if isinstance(rules, list):
            check_templates = [
                isinstance(c, Classes.ReactionTemplate) for c in rules
            ]
            assert all(
                check_templates
            ), """SubstructureNetwork.remove_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""
        else:
            assert isinstance(
                rules, list
            ), """SubstructureNetwork.remove_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""

        for temp in rules:
            self.remove_reaction_rule(temp)

    def add_reaction(self, reaction):
        """
        Add a reaction to the network.

        Parameters
        ----------
        reaction: NorthNet.Classes.Reaction

        Returns
        -------
        None
        """

        assert isinstance(
            reaction, Classes.Reaction
        ), """SubstructureNetwork.add_reaction:
            reaction arg should be a NorthNet Reaction object."""

        r_key = reaction.ReactionTemplate.ReactionSMARTS
        self.ReactionRules[r_key] = reaction

        for r_subst in reaction.ReactionTemplate.ReactantSubstructures:
            # connect to reaction
            if r_subst in self.Substructures:
                pass
            else:
                self.Substructures[r_subst] = Classes.Substructure(r_subst)

            self.Substructures[r_subst].ReactionParticipations.append(r_key)

            working_substruct = self.Substructures[r_subst].Mol

            # connect to compound
            for reac in reaction.Reactants:

                if reac in self.Compounds:
                    pass
                else:
                    self.Compounds[reac] = Classes.Compound(reac)

                compound = self.Compounds[reac].Mol

                if compound.HasSubstructMatch(working_substruct):
                    compound.ReactiveSubstructures.append(r_subst)
                    self.Substructures[r_subst].MatchingCompounds.append(reac)
                else:
                    pass

        for p_substruct in reaction.ReactionTemplate.ProductSubstructures:

            # connect to reaction
            if p_substruct in self.Substructures:
                pass
            else:
                self.Substructures[p_substruct] = Classes.Substructure(p_substruct)

            self.Substructures[p_substruct].ProducingReactions.append(r_key)

            working_substruct = self.Substructures[p_substruct].Mol

            # connect to compound
            for prod in reaction.Products:

                if prod in self.Compounds:
                    pass
                else:
                    self.Compounds[prod] = Classes.Compound(prod)

                compound = self.Compounds[prod].Mol

                if compound.HasSubstructMatch(working_substruct):
                    self.Compounds[prod].ReactantSubstructures.append(p_substruct)
                    self.Substructures[p_substruct].MatchingCompounds.append(prod)
                else:
                    pass

    def add_reactions(self, reactions):
        """
        Add reactions to the network.

        Parameters
        ----------
        reactions: list[NorthNet.Classes.Reaction]

        Returns
        -------
        None
        """

        if isinstance(reactions, list):
            check_reactions = [isinstance(c, Classes.Reaction) for c in reactions]
            assert all(
                check_reactions
            ), """SubstructureNetwork.add_reactions():
                reactions arg should be a list of NorthNet Reaction objects"""
        else:
            assert isinstance(
                reactions, list
            ), """SubstructureNetwork.add_reactions():
                reactions arg should be a list of NorthNet Reaction objects"""

        for r in reactions:
            if r.ReactionTemplate is None:
                pass
            elif r.ReactionTemplate.ReactionSMARTS in self.ReactionRules:
                pass
            else:
                self.add_reaction(r)

    def convert_to_networkx(self):
        """
        Converts NorthNet network object to networkx object.

        Parameters
        ----------

        Returns
        -------
        G: networkx.DiGraph
            Networkx version of the NorthNet SNetwork.
        """

        G = nx.DiGraph()

        # create some aliases for the substructures
        # (they cannot be SMARTS strings) used as node names
        substructure_aliases = {s: c for c, s in enumerate(self.Substructures)}

        for c in self.Compounds:
            compound_alias = self.Compounds[c].SMILES

            for i in self.Compounds[c].ReactiveSubstructures:
                G.add_edge(substructure_aliases[i], compound_alias)

        for r in self.ReactionRules:

            template = self.ReactionRules[r].ReactionTemplate

            transform_alias = template.Name

            for reac in template.ReactantSubstructures:
                G.add_edge(substructure_aliases[reac], transform_alias)

            for prod in template.ProductSubstructures:
                G.add_edge(transform_alias, substructure_aliases[prod])

        return G
