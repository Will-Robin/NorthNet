import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from NorthNet import Classes

class Substructure:
    '''
    Class to store substructures.
    '''
    def __init__(self,SMARTS):
        '''
        Parameters
        ----------
        SMARTS: str
            SMARTS corresponding to substructure.
        '''

        if isinstance(SMARTS, str):
            sys.exit('class Substructure: argument should be a SMARTS string.')

        self.Mol = Chem.MolFromSmarts(SMARTS)

        if Chem.MolToSmarts(self.Mol) == None:
            self.SMARTS = SMARTS
        else:
            self.SMARTS = Chem.MolToSmarts(self.Mol)

        self.In = []
        self.Out = []

class Compound:
    '''
    Class to store compound information.
    '''
    def __init__(self,SMILES):
        '''
        Parameters
        ----------
        mol: str
            SMILES corresponding to compound.
        '''

        if not isinstance(SMILES, str):
            sys.exit('class Compound: argument should be a SMILES string.')

        self.Mol = Chem.MolFromSmiles(SMILES)
        if self.Mol == None:
            self.SMILES = SMILES
        else:
            self.SMILES = Chem.MolToSmiles(self.Mol)

        self.In = []
        self.Out = []

        self.ReactiveSubstructures = None

class ReactionTemplate:
    '''
    Class for reaction templates.
    '''
    def __init__(self, name, reaction_SMARTS,
                        reactant_substructs, product_substructs):
        '''
        Parameters
        ----------
        name: str
            name for reaction.
        reaction_SMARTS: str
            reaction SMARTS string.

        reactant_substructs: list of strings
            reacting substructures

        product_substructs: list of strings
            substructures to which reactant substructures are converted.
        '''

        if not isinstance(name, str):
            sys.exit('''class ReactionTemplate:
        the name arg should be a name for the ReactionTemplate as a string.''')
        if not isinstance(reaction_SMARTS, str):
            sys.exit('''class ReactionTemplate:
            the reaction_SMARTS arg should be a reaction SMARTS string.''')
        if not isinstance(reactant_substructs, list):
            sys.exit('''class ReactionTemplate:
                    the reactant_substructs arg should be a SMARTS string.''')
        if not isinstance(product_substructs, list):
            sys.exit('''class ReactionTemplate:
                    the product_substructs arg should be a SMARTS string.''')

        self.Name = name
        if reaction_SMARTS != '':
            self.Reaction = AllChem.ReactionFromSmarts(reaction_SMARTS)
            self.ReactionSMARTS = AllChem.ReactionToSmarts(self.Reaction)
        else:
            self.Reaction = None
            self.ReactionSMARTS = ''

        self.ReactantSubstructures = reactant_substructs
        self.ProductSubstructures  = product_substructs

class Reaction:
    '''
    A class representing chemical reactions
    '''
    def __init__(self, rdkit_reaction, reaction_template = None, info = {}):
        '''

        Parameters
        ----------

        rdkit_reaction: rdkit Reaction object
            Reaction object for the reaction.
        reaction_template: NorthNet Reaction_Template object
            Reaction template which created the reaction.
        info: dict of dicts
            Dictionaries of dictionaries of information (e.g. database entries)

        Attributes
        ----------
        self.Reaction: rdkit reaction object
            Reaction in RDKit format
        self.ReactionSMILES: str
            Reaction SMILES string
        self.Reaction_Template: NorthNet ReactionTemplate or None
            Reaction template for the reaction. Defaults to None
        self.Data: dict of dicts (default to empty dict)
            Dictionaries of dictionaries of information (e.g. database entries)
        self.Reactants: list of str
            Reactant SMILES extracted from rdkit_reaction
        self.Products: list of str
            Product SMILES extracted from rdkit_reaction
        '''

        self.Reaction = rdkit_reaction
        self.ReactionSMILES = AllChem.ReactionToSmiles(rdkit_reaction)
        if reaction_template == None:
            self.ReactionTemplate = Classes.ReactionTemplate('none', '', [], [])
        else:
            self.ReactionTemplate = reaction_template
        self.Data = info

        self.Reactants = [Chem.MolToSmiles(x, canonical = True)
                                        for x in rdkit_reaction.GetReactants()]
        self.Products  = [Chem.MolToSmiles(x, canonical = True)
                                        for x in rdkit_reaction.GetProducts()]

class NetworkInput:
    '''
    Class to store inputs into reaction network
    '''
    def __init__(self,id):
        '''
        Designed to behave like a compound object
        self.SMILES does not actually hold a valid SMILES string,
        it stores the input ID code.

        Parameters
        ----------
        id: str
            token for reaction input
            should follow the convention
            SMILES_#0
        '''

        if not isinstance(id, str):
            sys.exit('''class NetworkInput:
                        id arg must be str like SMILES_#0.''')

        self.Mol = Chem.MolFromSmiles(id.split('_')[0])

        # Note that the self.SMILES does not
        # actually hold a valid SMILES string,
        # it stores the input ID code
        self.SMILES = id

        self.In = None
        self.Out = []

        self.ReactiveSubstructures = None

class NetworkOutput:
    '''
    Class to store inputs into reaction network
    '''
    def __init__(self,id):
        '''
        Designed to behave like a compound object
        self.SMILES does not actually hold a valid SMILES string,
        it stores the input ID code.

        Parameters
        ----------
        id: str
            token for reaction input
            should follow the convention
            SMILES_#0
        '''

        if not isinstance(id, str):
            sys.exit('''class NetworkOutput:
                    id arg must be str like SMILES_#0.''')

        self.Mol = None

        # Note that the self.SMILES does not
        # actually hold a valid SMILES string,
        # it stores the input ID code
        self.SMILES = id

        self.In = []
        self.Out = None

        self.ReactiveSubstructures = None

class ReactionInput:
    '''
    Class to store inputs as reactions into reaction network
    '''
    def __init__(self, reaction_input_string):
        '''
        Designed to behave like a Reaction object
        self.SMILES does not actually hold a valid SMILES string,

        Parameters
        ----------
        reaction_input_string: str
            token for reaction input
            should follow the convention
            SMILES_#0>>SMILES
        '''

        if not isinstance(id, str):
            sys.exit('''class ReactionInput:
                    reaction_input_string arg must be str like SMILES_#0>>SMILES.''')

        self.Reaction = None
        self.ReactionSMILES = reaction_input_string
        self.Reactants, self.Products = self.reactants_products_from_string(
                                                        reaction_input_string)

        self.CompoundInput = self.Products[0]
        self.InputID = self.Reactants[0]

        self.ReactionTemplate = None
        self.Data = {}

    def reactants_products_from_string(self, reaction_smiles):
        split_rxn_smiles = reaction_smiles.split('>>')
        reactants = [x for x in split_rxn_smiles[0].split('.') if x != '']
        products = [x for x in split_rxn_smiles[1].split('.') if x != '']
        return reactants, products

class ReactionOutput:
    '''
    Class to store outputs as reactions into reaction network
    '''
    def __init__(self, reaction_output_string):
        '''
        Designed to behave like a Reaction object
        self.SMILES does not actually hold a valid SMILES string,

        Parameters
        ----------
        reaction_input_string: str
            token for reaction input
            should follow the convention
            SMILES_#0>>SMILES
        '''
        if not isinstance(reaction_output_string, str):
            sys.exit('''class ReactionOutput:
                    reaction_output_string must be string like SMILES>>#0''')
        self.Reaction = None
        self.ReactionSMILES = reaction_output_string
        self.Reactants, self.Products = self.reactants_products_from_string(
                                                        reaction_output_string)

        self.OutputID = self.Products[0]
        self.CompoundOutput = self.Reactants[0]

        self.ReactionTemplate = None
        self.Data = {}

    def reactants_products_from_string(self, reaction_smiles):
        split_rxn_smiles = reaction_smiles.split('>>')
        reactants = [x for x in split_rxn_smiles[0].split('.') if x != '']
        products = [x for x in split_rxn_smiles[1].split('.') if x != '']
        return reactants, products

class Network:
    def __init__(self, reactions, name, description):
        '''
        A network objects consisting of Edges (reactions) and nodes (species).
        The Network object is initialised with a list of Reaction objects. If
        the list is empty, then the network is initialised as an empty network.
        '''

        self.Name = name
        self.Description = description

        '''1. create a dictionary of lists of Reaction objects, which are
        organised by their reaction SMILES as keys. '''
        self.NetworkReactions = {}

        '''2. Create dictionary of Compound objects used as reactants
        and products from the reaction objects (SMILES used as keys).'''
        self.NetworkCompounds = {}

        '''3. Create dictionary to store network inputs'''
        self.NetworkInputs = {}

        '''4. Create dictionary of network outputs'''
        self.NetworkOutputs = {}

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_compounds(self,compounds):
        '''
        Add list of NorthNet Compound objects to the Network

        Parameters
        ----------
        compounds: list of NorthNet Compound objects
            compounds to be added
        '''
        for n in compounds:
            if n.SMILES in self.NetworkCompounds:
                pass
            else:
                self.NetworkCompounds[n.SMILES] = n

    def remove_compounds(self, compounds):
        '''
        Remove list of compounds and the reactions in which they are involved
        from the network.

        Parameters
        ----------
        compounds: list of NorthNet Compound objects
            compounds to be removed
        '''

        remove_reactions = []
        for c in compounds:
            remove_reactions.extend(self.NetworkCompounds[c.SMILES].In)
            remove_reactions.extend(self.NetworkCompounds[c.SMILES].Out)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions([self.NetworkReactions[r]
                                                    for r in remove_reactions])

        for c in compounds:
            del self.NetworkCompounds[c.SMILES]

    def add_reaction(self, reaction):
        '''
        Adds a reaction and its associated reactants and products into the
        Network.

        Parameters
        ----------
        reaction: NorthNet Reaction object
            reaction to be added
        '''

        # check if the reaction is in NetworkReactions (avoid overwriting)
        # shouldn't be an issue if all reaction data is encapsulated properly
        # i.e. one reaction SMILES string to one set of data
        if reaction.ReactionSMILES in self.NetworkReactions:
            pass
        else:
            # add the reaction into NetworkReactions
            self.NetworkReactions[reaction.ReactionSMILES] = reaction

            for a in reaction.Reactants:
                if a in self.NetworkCompounds:
                    # if reaction is already connected to the reactant,
                    # pass. This should not happen if all reaction data
                    # is encapsulated properly
                    if reaction.ReactionSMILES in self.NetworkCompounds[a].Out:
                        pass
                    # connect the reactant to the reaction
                    else:
                        self.NetworkCompounds[a].Out.append(
                                                        reaction.ReactionSMILES)
                # add the reactant into NetworkCompounds
                # and connect the reactant to the reaction
                else:
                    self.NetworkCompounds[a] = Classes.Compound(a)
                    self.NetworkCompounds[a].Out.append(reaction.ReactionSMILES)

            for b in reaction.Products:
                if b in self.NetworkCompounds:
                    # if reaction is already connected to the product,
                    # pass. This should not happen if all reaction data
                    # is encapsulated properly
                    if reaction.ReactionSMILES in self.NetworkCompounds[b].In:
                        pass
                    # connect the reaction to the product
                    else:
                        self.NetworkCompounds[b].In.append(
                                                        reaction.ReactionSMILES)
                # add the product into NetworkCompounds
                # and connect the reaction to the product
                else:
                    self.NetworkCompounds[b] = Classes.Compound(b)
                    self.NetworkCompounds[b].In.append(reaction.ReactionSMILES)

    def add_reactions(self,reactions):
        '''
        Use the standardised strings information in the Reaction objects
        to build them into the Network.

        Parameters
        ----------
        reactions: list of NorthNet Reaction objects
            reactions to be added to the Network
        '''
        for r in reactions:
            self.add_reaction(r)

    def remove_reactions(self, remove_reactions):
        '''
        Remove a list of reactions from the Network.

        Parameters
        ----------
        remove_reactions: list of NorthNet Reaction objects
            reactions to be removed.
        '''

        for r in remove_reactions:
            r_key = r.ReactionSMILES
            for rc in self.NetworkReactions[r_key].Reactants:
                self.NetworkCompounds[rc].Out.remove(r_key)
            for p in self.NetworkReactions[r_key].Products:
                self.NetworkCompounds[p].In.remove(r_key)

            del self.NetworkReactions[r_key]

    def add_input(self, input):
        '''
        Add a NetworkInput to the Network

        Parameters
        ----------
        input: NorthNet NetworkInput object
            Input to be added
        '''
        if input.InputID in self.NetworkInputs:
            pass
        else:
            self.NetworkReactions[input.ReactionSMILES] = input

            if input.CompoundInput in self.NetworkCompounds:
                pass
            else:
                self.NetworkCompounds[input.CompoundInput] = input

            self.NetworkCompounds[input.CompoundInput].In.append(
                                                        input.ReactionSMILES)

            if input.InputID in self.NetworkInputs:
                pass
            else:
                self.NetworkInputs[input.InputID] = input

            self.NetworkInputs[input.InputID].Out.append(input.ReactionSMILES)

    def add_inputs(self, inputs):
        '''
        For adding NetworkInput to the network

        Parameters
        ----------
        inputs: list of NortNet NetworkInput objects
            Inputs to be added to the network
        '''
        for i in inputs:
            self.add_input(i)

    def add_output(self, output):
        '''
        Add a NetworkOutput to the Network

        Parameters
        ----------
        output NetworkOutput object
            Output to be added
        '''
        if output.OutputID in self.NetworkOutputs:
            pass
        else:
            if output.CompoundOutput in self.NetworkCompounds:
                self.NetworkReactions[output.ReactionSMILES] = output
                self.NetworkCompounds[output.CompoundOutput].Out.append(
                                                        output.ReactionSMILES)
            else:
                # The compound is not currently in the network, so cannot
                # be an output
                pass

            if output.OutputID in self.NetworkOutputs:
                pass
            else:
                self.NetworkOutputs[output.OutputID] = output

            self.NetworkOutputs[output.OutputID].In.append(
                                                    output.ReactionSMILES)

    def add_outputs(self, outputs):
        '''
        For adding NetworkOutput to the network

        Parameters
        ----------
        outputs: list of NortNet NetworkOutput objects
            Outputs to be added to the network
        '''
        for o in outputs:
            self.add_output(o)

    def get_reaction(self, reaction):
        '''
        Convenience class for getting a reaction using a key

        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.NetworkReactions[reaction]

    def get_reactants(self, reaction):
        '''
        Conveniently get the reactants of a reaction
        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.get_reaction(reaction).Reactants

    def get_products(self, reaction):
        '''
        Conveniently get the products of a reaction
        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.get_reaction(reaction).Products

    def get_reaction_template(self, reaction):
        '''
        Conveniently get the ReactionTemplate of a reaction
        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.get_reaction(reaction).ReactionTemplate

    def get_reaction_SMARTS(self, reaction):
        '''
        Conveniently get the Reaction SMARTS of a reaction
        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.get_reaction_template(reaction).ReactionSMARTS

    def get_reaction_name(self, reaction):
        '''
        Conveniently get the Name of a reaction
        Parameters
        ----------
        reaction: str
            Key in self.NetworkReactions
        '''
        return self.get_reaction_template(reaction).Name

    def convert_to_networkx(self):
        '''
        Converts NorthNet network object to networkx object.

        Returns
        -------
        G: networkx DiGraph object
            Networkx version of the NorthNet network.
        '''
        import networkx as nx

        G = nx.DiGraph()

        for node in self.NetworkCompounds:
            G.add_node(node)

        for r in self.NetworkReactions:
            for sr in self.NetworkReactions[r].Reactants:
                for sp in self.NetworkReactions[r].Products:
                    G.add_edge(sr,r)
                    G.add_edge(r,sp)

        return G

class SubstructureNetwork:
    '''
    A network designed to show the relationship between functional
    group transformations.
    '''
    def __init__(self,reactions,name,description):
        '''
        reactions: list
            List of reaction objects with extracted functional group
            transformations.
        name: str
            A name for the network.
        description: str
            A description of the network.
        '''

        '''
        Basic idea:
            Nodes will be individual reaction centres and functional group
            transformations. Edges will connect the reaction centres to
            functional group transformations.
        '''

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
        reactions: NorthNet Reaction object
        '''
        # if reaction.ReactionTemplate != None:

        r_key = reaction.ReactionTemplate.ReactionSMARTS
        self.SNetworkTemplates[r_key] = reaction

        for r_subst in reaction.ReactionTemplate.ReactantSubstructures:
            # connect to reaction
            if r_subst in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[r_subst] = Classes.Substructure(r_subst)

            self.SNetworkSubstructs[r_subst].Out.append(r_key)

            working_substruct = self.SNetworkSubstructs[r_substruct].Mol
            # connect to compound
            for reac in reactant.Reactants:

                if reac in self.SNetworkCompounds:
                    pass
                else:
                    self.SNetworkCompounds[reac] = Classes.Compound(reac)

                compound = self.SNetworkCompounds[reac].Mol

                if compound.HasSubstructMatch(working_substruct):
                    compound.Out.append(r_substruct)
                    self.SNetworkSubstructs[r_substruct].In.append(reac)
                else:
                    pass

        for p_substruct in r.ReactionTemplate.ProductSubstructures:

            # connect to reaction
            if p_substruct in self.SNetworkSubstructs:
                pass
            else:
                self.SNetworkSubstructs[p_substruct] = Classes.Substructure(
                                                                    p_substruct)

            self.SNetworkSubstructs[p_substruct].In.append(r_key)

            working_substruct = self.SNetworkSubstructs[p_substruct].Mol
            # connect to compound
            for prod in r.Products:

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
        Parameters
        ----------
        reactions: list of NorthNet Reaction objects
        '''

        for r in reactions:
            if r.ReactionTemplate == None:
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

            for i in snetwork.SNetworkCompounds[n].In:
                G.add_edge(substructure_aliases[i], compound_alias)
            for o in snetwork.SNetworkCompounds[n].Out:
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
