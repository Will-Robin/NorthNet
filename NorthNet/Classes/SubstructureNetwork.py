import sys
from NorthNet import Classes

class SubstructureNetwork:
    '''
    TODO: fill out with similar methods as for the Network object.

    A network designed to show the relationship between functional
    group transformations.

    The network has three kinds of node: reactions, compounds, and
    substructures. Substructures connect reactions to compounds. Compounds and
    reactions can only link to substructures.
    '''
    def __init__(self,reactions,name,description):
        '''
        reactions: list
            List of NortNet Reaction objects with extracted functional group
            transformations.
        name: str
            A name for the network.
        description: str
            A description of the network.
        '''

        if isinstance(reactions, list):
            check_reactions = [isinstance(r, Classes.Reaction) for r in reactions]
            if not all(check_reactions):
                sys.exit('''class SubstructureNetwork: the reactions arg must be a list of
                NorthNet Reaction objects.''')
        else:
            sys.exit('''class SubstructureNetwork: the reactions arg must be a list of
            NorthNet Reaction objects.''')

        if not isinstance(name, str):
            sys.exit('''class SubstructureNetwork: name arg should be a string.''')
        if not isinstance(name, str):
            sys.exit('''class SubstructureNetwork: description arg should be a string.''')

        self.Name = name
        self.Description = description

        self.SNetworkSubstructs = {}

        self.SNetworkTemplates = {}

        self.SNetworkCompounds = {}

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_reaction(self, reaction):
        '''
        Add a reaction to the network.
            
        Parameters
        ----------
        reactions: NorthNet Reaction object
        '''
        if not isinstance(reaction, Classes.Reaction):
            sys.exit('''SubstructureNetwork.add_reaction: reaction arg should be a NorthNet
            Reaction object.''')

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
                self.SNetworkSubstructs[p_substruct] = Classes.Substructure(
                                                                    p_substruct)

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
        '''
        Add reactions to the network.

        Parameters
        ----------
        reactions: list of NorthNet Reaction objects
        '''
        
        if isinstance(reactions, list):
            check_reactions = [isinstance(c, Classes.Reaction) for c in reactions]
            if not all(check_reactions):
                sys.exit('''SubstructureNetwork.add_reactions():
                reactions arg should be a list of NorthNet Reaction objects''')
        else:
            sys.exit('''SubstructureNetwork.add_reactions():
            reactions arg should be a list of NorthNet Reaction objects''')

        for r in reactions:
            if r.ReactionTemplate is None:
                pass
            elif r.ReactionTemplate.ReactionSMARTS in self.SNetworkTemplates:
                pass
            else:
                self.add_reaction(r)

    def convert_to_networkx(self):
        '''
        Converts NorthNet network object to networkx object.

        Parameters
        ----------
        Returns
        -------
        G: networkx DiGraph object
            Networkx version of the NorthNet SNetwork.
        '''

        import networkx as nx

        G = nx.DiGraph()
        # create some aliases for the substructures
        # (they cannot be SMARTS strings) used as node names
        substructure_aliases = {s:c
                                for c,s in enumerate(self.SNetworkSubstructs)}

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
