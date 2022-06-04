from NorthNet import Classes


class SubstructureNetwork:
    """
    TODO: fill out with similar methods as for the Network object.

    A network designed to show the relationship between functional
    group transformations.

    The network has three kinds of node: reactions, compounds, and
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
        SNetworkSubstructs: dict
            Dictionary containing NorthNet.Classes.Substructure objects keyed
            by SMARTS.
        SNetworkTemplates: dict
            Dictionary containing NorthNet.Classes.ReactionTemplate objects keyed
            by reaction SMARTS.
        SNetworkCompounds: dict
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
        self.SNetworkSubstructs = dict()
        self.SNetworkTemplates = dict()
        self.SNetworkCompounds = dict()

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_compound(self, compound):
        """
        Add a compound to the SubstructureNetwork.

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

        if compound.SMILES not in self.SNetworkCompounds:
            self.SNetworkCompounds[compound.SMILES] = compound

    def add_compounds(self, compounds):
        """
        Add list of NorthNet Compound objects to the SubstructureNetwork

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

    def remove_compounds(self, compounds):
        """
        Remove list of compounds from the SubstructureNetwork.

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

        check_substructs = []
        for compound in compounds:
            for substruct in self.SNetworkCompounds[compound].In:
                self.SNetworkSubstructs[substruct].Out.remove(compound)
                check_substructs.append(substruct)
            for substruct in self.SNetworkCompounds[compound].Out:
                self.SNetworkSubstructs[substruct].In.remove(compound)
                check_substructs.append(substruct)

        for substruct in check_substructs:

            substr_obj = self.SNetworkSubstructs[substruct]

            if len(substr_obj.Out) == 0 and len(substr_obj.In) == 0:
                remove_templates = []

                for r_t in substr_obj.In:
                    remove_templates.append(self.SNetworkTemplates[r_t])
                for r_t in substr_obj.Out:
                    remove_templates.append(self.SNetworkTemplates[r_t])

                remove_templates = list(set(remove_templates))
                templ_removals = [self.SNetworkTemplates[t] for t in remove_templates]

                self.remove_substructure(substr_obj)
                self.remove_templates(templ_removals)

        for compound in compounds:
            del self.SNetworkCompounds[compound.SMILES]

    def remove_substructure(self, substructure):
        """
        Remove a substructure from the SubstructureNetwork.

        Parameters
        ----------
        substructure: NorthNet.Classes.Substructure

        Returns
        -------
        None
        """

        assert isinstance(
            substructure, Classes.Substructure
        ), """SubstructureNetwork.remove_substructure:
            substructure arg should be a NorthNet Substructure object."""

        del self.SNetworkSubstructs[substructure.SMARTS]

    def remove_substructures(self, substructures):
        """
        Remove a list of substructures from the SubstructureNetwork.

        Parameters
        ----------
        substructures: list[NorthNet.Classes.Substructure]

        Returns
        -------
        None
        """

        if isinstance(substructures, list):
            check_subtructures = [
                isinstance(c, Classes.Substructure) for c in substructures
            ]
            assert all(
                check_subtructures
            ), """SubstructureNetwork.remove_substructures():
                substructures arg should be a list of NorthNet Substructure objects"""
        else:
            assert isinstance(
                substructures, list
            ), """SubstructureNetwork.remove_substructures():
                substructures arg should be a list of NorthNet Substructure objects"""

        for substruct in substructures:
            self.remove_substructure(substruct)

    def add_template(self, template):
        """
        Add a reaction template to the network.

        Parameters
        ----------
        template: NorthNet.Classes.ReactionTemplate

        Returns
        -------
        None
        """

        assert isinstance(
            template, Classes.ReactionTemplate
        ), """SubstructureNetwork.add_template:
            template arg should be a NorthNet ReactionTemplate object."""

        r_key = template.ReactionSMARTS

        for r_subst in template.ReactantSubstructures:
            # connect to reaction
            if r_subst in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[r_subst] = Classes.Substructure(r_subst)

            self.SNetworkSubstructs[r_subst].Out.append(r_key)

            working_substruct = self.SNetworkSubstructs[r_subst].Mol

            # connect to compounds
            for comp in self.SNetworkCompounds:
                compound = self.SNetworkCompounds[comp].Mol
                if compound.HasSubstructMatch(working_substruct):
                    self.SNetworkCompounds[comp].In.append(r_subst)
                    self.SNetworkSubstructs[r_subst].Out.append(comp)

        for p_substruct in template.ProductSubstructures:
            # connect to reaction
            if p_substruct in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[p_substruct] = Classes.Substructure(p_substruct)

            self.SNetworkSubstructs[p_substruct].In.append(r_key)

            working_substruct = self.SNetworkSubstructs[p_substruct].Mol

            # connect to compounds
            for comp in self.SNetworkCompounds:
                compound = self.SNetworkCompounds[comp].Mol
                if compound.HasSubstructMatch(working_substruct):
                    self.SNetworkCompounds[comp].Out.append(p_substruct)
                    self.SNetworkSubstructs[p_substruct].In.append(comp)

    def add_templates(self, templates):
        """
        Add reaction templates to the network.

        Parameters
        ----------
        reactions: list[NorthNet.Classes.ReactionTemplate]

        Returns
        -------
        None
        """

        if isinstance(templates, list):
            check_templates = [
                isinstance(c, Classes.ReactionTemplate) for c in templates
            ]
            assert all(
                check_templates
            ), """SubstructureNetwork.add_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""
        else:
            assert isinstance(
                templates, list
            ), """SubstructureNetwork.add_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""

        for t in templates:
            if t.ReactionSMARTS in self.SNetworkTemplates:
                pass
            else:
                self.add_template(t)

    def remove_template(self, template):
        """
        Remove a reaction template from a SubstructureNetwork.

        Parameters
        ----------
        template: NorthNet.Classes.ReactionTemplate

        Returns
        -------
        None
        """

        assert isinstance(
            template, Classes.ReactionTemplate
        ), """SubstructureNetwork.remove_template:
            template arg should be a NorthNet ReactionTemplate object."""

        del self.SNetworkTemplates[template.ReactionSMARTS]

    def remove_templates(self, templates):
        """
        Remove a list of ReactionTemplates from the SubstructureNetwork.

        Parameters
        ----------
        templates: list[NorthNet.Classes.ReactionTemplate]

        Returns
        -------
        None
        """

        if isinstance(templates, list):
            check_templates = [
                isinstance(c, Classes.ReactionTemplate) for c in templates
            ]
            assert all(
                check_templates
            ), """SubstructureNetwork.remove_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""
        else:
            assert isinstance(
                templates, list
            ), """SubstructureNetwork.remove_templates():
                templates arg should be a list of NorthNet ReactionTemplate objects"""

        for temp in templates:
            self.remove_template(temp)

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
        self.SNetworkTemplates[r_key] = reaction

        for r_subst in reaction.ReactionTemplate.ReactantSubstructures:
            # connect to reaction
            if r_subst in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[r_subst] = Classes.Substructure(r_subst)

            self.SNetworkSubstructs[r_subst].Out.append(r_key)

            working_substruct = self.SNetworkSubstructs[r_subst].Mol

            # connect to compound
            for reac in reaction.Reactants:

                if reac in self.SNetworkCompounds:
                    pass
                else:
                    self.SNetworkCompounds[reac] = Classes.Compound(reac)

                compound = self.SNetworkCompounds[reac].Mol

                if compound.HasSubstructMatch(working_substruct):
                    compound.Out.append(r_subst)
                    self.SNetworkSubstructs[r_subst].In.append(reac)
                else:
                    pass

        for p_substruct in reaction.ReactionTemplate.ProductSubstructures:

            # connect to reaction
            if p_substruct in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[p_substruct] = Classes.Substructure(p_substruct)

            self.SNetworkSubstructs[p_substruct].In.append(r_key)

            working_substruct = self.SNetworkSubstructs[p_substruct].Mol
            # connect to compound
            for prod in reaction.Products:

                if prod in self.SNetworkCompounds:
                    pass
                else:
                    self.SNetworkCompounds[prod] = Classes.Compound(prod)

                compound = self.SNetworkCompounds[prod].Mol

                if compound.HasSubstructMatch(working_substruct):
                    self.SNetworkCompounds[prod].In.append(p_substruct)
                    self.SNetworkSubstructs[p_substruct].Out.append(prod)
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
            elif r.ReactionTemplate.ReactionSMARTS in self.SNetworkTemplates:
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

        import networkx as nx

        G = nx.DiGraph()
        # create some aliases for the substructures
        # (they cannot be SMARTS strings) used as node names
        substructure_aliases = {s: c for c, s in enumerate(self.SNetworkSubstructs)}

        for s in self.SNetworkSubstructs:
            G.add_node(substructure_aliases[s])

        for c in self.SNetworkCompounds:
            compound_alias = self.SNetworkCompounds[c].SMILES

            G.add_node(compound_alias)

            for i in self.SNetworkCompounds[c].In:
                G.add_edge(substructure_aliases[i], compound_alias)

            for o in self.SNetworkCompounds[c].Out:
                G.add_edge(compound_alias, substructure_aliases[o])

        for r in self.SNetworkTemplates:
            template = self.SNetworkTemplates[r].ReactionTemplate
            transform_alias = template.Name

            G.add_node(transform_alias)

            for reac in template.ReactantSubstructures:
                G.add_edge(substructure_aliases[reac], transform_alias)

            for prod in template.ProductSubstructures:
                G.add_edge(transform_alias, substructure_aliases[prod])

        return G
