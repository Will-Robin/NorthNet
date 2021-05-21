from rdkit.Chem import AllChem
from rdkit import Chem
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
        self.Mol = Chem.MolFromSmiles(SMILES)
        if self.Mol == None:
            self.SMILES = SMILES
        else:
            self.SMILES = Chem.MolToSmiles(self.Mol)

        self.In = []
        self.Out = []
        self.Block = False
        self.EdgeLog = []

        self.ReactiveSubstructures = None

class Reaction_Template:
    '''
    Class for reaction templates.
    '''
    def __init__(self,name,reaction_SMARTS, reactant_substructs, product_substructs):
        '''
        Parameters
        ----------
        name: str
            name for reaction.
        reaction_SMARTS: list
            List of reaction SMARTS strings.
        '''

        self.Name = name
        self.Reaction = AllChem.ReactionFromSmarts(reaction_SMARTS)
        self.ReactionSMARTS = AllChem.ReactionToSmarts(self.Reaction)
        self.ReactantSubstructures = reactant_substructs
        self.ProductSubstructures  = product_substructs

class Generated_Reaction:
    '''
    An object designed to store generated reactions and their characteristics.
    '''
    def __init__(self,rdkit_reaction,reaction_template):
        '''
        Parameters
        ----------
        rdkit_reaction: rdkit Reaction object
            Reaction object for the reaction.
        reaction_template: NorthNet Reaction_Template object
            Reaction template which created the reaction.
        '''

        self.Reaction = rdkit_reaction
        self.ReactionSMILES = AllChem.ReactionToSmiles(rdkit_reaction)
        self.Reaction_Template = reaction_template

class Reaction_Database_Entry:
        '''
        Class for storing information on a single reaction from a database.
        '''
        def __init__(self, header, line):
            '''
            Parameters
            ----------
            header: list
                Header for field entries.
            line: list
                List of field entries from database.
            '''

            self.Info = {k:v for k,v in zip(header,line)}
            self.ReactionSMILES = line[header.index('Reaction')].strip(',')

        def recompile(self):
            '''Creates a modified database entry by modifiying the original with the
            self.Info.'''

            newline = []
            for k in [*self.Info]:
                newline.append(self.Info[k])

            return newline

        def get_header(self):
            '''
            For getting the field entries for the database entry.
            '''

            newline = []
            for k in [*self.Info]:
                newline.append(k)

            return newline

class Reaction:
    '''
    A class which unites Generated_Reaction and Reaction_Database_Entry classes
    with the same reaction SMILES so they are combined in the view of a Network.

    Could use multiple inheritance, but this seems simpler.
    '''
    def __init__(self, reaction_object):
        '''
        Parameters
        ----------
        reaction_object: NorthNet Generated_Reaction or Reaction_Database_Entry objects
            Prototypical reaction to initialise the object.

        Attributes
        ----------
        self.Database_Entries: list
            List of Reaction_Database_Entry objects
        self.Generation_Details: list
            list of Generated_Reaction objects
        '''

        if reaction_object.__class__.__name__ == "Reaction_Database_Entry":
            self.Database_Entries = [reaction_object]
            self.Generation_Details = []

            spl = reaction_object.ReactionSMILES.split('>>')

            self.Reactants = [r for r in spl[0].split('.') if Chem.MolFromSmiles(r) != None and r != ""]
            self.Products  = [p for p in spl[1].split('.') if Chem.MolFromSmiles(p) != None and p != ""]

            self.ReactionTemplate = None
            self.MappedReaction = None

            self.update_reaction()

        if reaction_object.__class__.__name__ == "Generated_Reaction":
            self.Database_Entries = []
            self.Generation_Details = [reaction_object]
            self.Reactants = [Chem.MolToSmiles(r) for r in reaction_object.Reaction.GetReactants()]
            self.Products  = [Chem.MolToSmiles(p) for p in reaction_object.Reaction.GetProducts()]

            self.ReactionTemplate = reaction_object.Reaction_Template
            self.MappedReaction = None # If mapping can be introduced inot the reaction generation process. If it is, this can be changed accordingly.
            self.update_reaction()

        if reaction_object.__class__.__name__ == "str":
            self.Database_Entries = []
            self.Generation_Details = []
            self.Reactants, self.Products = self.reactants_products_from_string(reaction_object)
            self.ReactionTemplate = None
            self.MappedReaction = None # If mapping can be introduced inot the reaction generation process. If it is, this can be changed accordingly.
            self.update_reaction()

    def add_reaction_entry(self, new_entry):
        if new_entry.__class__.__name__ == "Reaction_Database_Entry":
            self.Database_Entries.append(new_entry)
        if new_entry.__class__.__name__ == "Generated_Reaction":
            self.Generation_Details.append(new_entry)

    def reactants_products_from_string(self, reaction_smiles):
        split_rxn_smiles = reaction_smiles.split('>>')
        reactants = split_rxn_smiles[0].split('.')
        products = split_rxn_smiles[1].split('.')
        return reactants, products

    def update_reaction(self):
        '''
        Updates ReactionSMILES and Species attributes according to self.Reactants
        and self.Products.
        '''
        rxn = AllChem.ChemicalReaction() # Create an empty chemical reaction
        [rxn.AddReactantTemplate(Chem.MolFromSmiles(r)) for r in self.Reactants]
        [rxn.AddProductTemplate(Chem.MolFromSmiles(p)) for p in self.Products]
        self.Reaction = rxn
        self.ReactionSMILES = AllChem.ReactionToSmiles(rxn)

        for d in self.Database_Entries:
            d.ReactionSMILES = self.ReactionSMILES
            d.Info["Reaction"] = self.ReactionSMILES

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

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_reactions(self,reactions):
        for r in reactions:

            if r.ReactionSMILES in self.NetworkReactions:
                self.NetworkReactions[r.ReactionSMILES].add_reaction_entry(r)

            else:
                self.NetworkReactions[r.ReactionSMILES] = r

                for a in r.Reactants:
                    if a in self.NetworkCompounds:
                        if r.ReactionSMILES in self.NetworkCompounds[a].Out:
                            pass
                        else:
                            self.NetworkCompounds[a].Out.append(r.ReactionSMILES)
                    else:
                        self.NetworkCompounds[a] = Classes.Compound(a)
                        self.NetworkCompounds[a].Out = [r.ReactionSMILES]

                for b in r.Products:
                    if b in self.NetworkCompounds:
                        if r.ReactionSMILES in self.NetworkCompounds[b].In:
                            pass
                        else:
                            self.NetworkCompounds[b].In.append(r.ReactionSMILES)
                    else:
                        self.NetworkCompounds[b] = Classes.Compound(b)
                        self.NetworkCompounds[b].In = [r.ReactionSMILES]

    def add_compounds(self,compounds):
        for n in compounds:
            if n.SMILES in self.NetworkCompounds:
                pass
            else:
                self.NetworkCompounds[n.SMILES] = n

    def remove_compounds(self, compounds):

        remove_reactions = []
        for c in compounds:
            remove_reactions.extend(self.NetworkCompounds[c].In)
            remove_reactions.extend(self.NetworkCompounds[c].Out)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions(remove_reactions)

        for c in compounds:
            del self.NetworkCompounds[c]

    def remove_reactions(self, remove_reactions):
        for r in remove_reactions:
            for rc in self.NetworkReactions[r].Reactants:
                self.NetworkCompounds[rc].Out.remove(r)
            for p in self.NetworkReactions[r].Products:
                self.NetworkCompounds[p].In.remove(r)

            del self.NetworkReactions[r]


class Substructure_Network:
    '''
    A network designed to show the relationship between functional
    group transformations.
    '''
    def __init__(self,reactions,name,description):
        '''
        reactions: list
            List of reaction objects with extracted functional group transformations.
        name: str
            A name for the network.
        description: str
            A description of the network.
        '''

        '''
        Basic idea:
            Nodes will be individual reaction centres and functional group transformations
            Edges will connect the reaction centres to functional group transformations.
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

    def add_reactions(self, reactions):

        for r in reactions:

            r_key = r.ReactionTemplate.ReactionSMARTS

            if r_key in self.SNetworkTemplates:
                pass
            else:
                self.SNetworkTemplates[r_key] = r

            for r_substruct in r.ReactionTemplate.ReactantSubstructures:
                # connect to reaction
                if r_substruct in self.SNetworkSubstructs:
                    self.SNetworkSubstructs[r_substruct].Out.append(r_key)
                else:
                    self.SNetworkSubstructs[r_substruct] = Classes.Substructure(r_substruct)
                    self.SNetworkSubstructs[r_substruct].Out.append(r_key)

                # connect to compound
                for reac in r.Reactants:
                    if reac in self.SNetworkCompounds:
                        if self.SNetworkCompounds[reac].Mol.HasSubstructMatch(self.SNetworkSubstructs[r_substruct].Mol):
                            self.SNetworkCompounds[reac].Out.append(r_substruct)
                            self.SNetworkSubstructs[r_substruct].In.append(reac)
                        else:
                            pass
                    else:
                        self.SNetworkCompounds[reac] = Classes.Compound(reac)

                        if self.SNetworkCompounds[reac].Mol.HasSubstructMatch(self.SNetworkSubstructs[r_substruct].Mol):
                            self.SNetworkCompounds[reac].Out.append(r_substruct)
                            self.SNetworkSubstructs[r_substruct].In.append(reac)
                        else:
                            pass

            for p_substruct in r.ReactionTemplate.ProductSubstructures:
                # connect to reaction
                if p_substruct in self.SNetworkSubstructs:
                    self.SNetworkSubstructs[p_substruct].In.append(r_key)
                else:
                    self.SNetworkSubstructs[p_substruct] = Classes.Substructure(p_substruct)
                    self.SNetworkSubstructs[p_substruct].In.append(r_key)

                # connect to compound
                for prod in r.Products:
                    if p_substruct == "[Ch:3][O:4]" and prod in self.SNetworkCompounds:
                        pass
                        #print(prod,p_substruct,self.SNetworkCompounds[prod].Mol.HasSubstructMatch(self.SNetworkSubstructs[p_substruct].Mol))
                    if prod in self.SNetworkCompounds:
                        if self.SNetworkCompounds[prod].Mol.HasSubstructMatch(self.SNetworkSubstructs[p_substruct].Mol):
                            self.SNetworkCompounds[prod].In.append(p_substruct)
                            self.SNetworkSubstructs[p_substruct].Out.append(prod)
                        else:
                            pass
                    else:
                        self.SNetworkCompounds[prod] = Classes.Compound(prod)
                        if self.SNetworkCompounds[prod].Mol.HasSubstructMatch(self.SNetworkSubstructs[p_substruct].Mol):
                            self.SNetworkCompounds[prod].In.append(p_substruct)
                            self.SNetworkSubstructs[p_substruct].Out.append(prod)
                        else:
                            pass
