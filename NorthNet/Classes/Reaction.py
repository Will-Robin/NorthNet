from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

from NorthNet import Classes

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
        info: dict
            Dictionary of information (e.g. database entries)

        Attributes
        ----------
        self.Reaction: rdkit reaction object
            Reaction in RDKit format
        self.ReactionSMILES: str
            Reaction SMILES string
        self.Reaction_Template: NorthNet ReactionTemplate or None
            Reaction template for the reaction. Defaults to None
        self.info: dict
            Dictionary of information (e.g. database entries)
        '''

        assert isinstance(rdkit_reaction, rdChemReactions.ChemicalReaction), \
            '''class Reaction:
            rdkit_reaction arg should be an RDKit Reaction object.'''

        if reaction_template is not None:
            assert isinstance(reaction_template, Classes.ReactionTemplate), \
                '''class Reaction:
                reaction_template should be None or a NortNet
                ReactionTemplate object.'''

        assert isinstance(info, dict), \
            '''class Reaction: info arg should be a dict.'''

        self.Reaction = rdkit_reaction
        self.ReactionSMILES = AllChem.ReactionToSmiles(rdkit_reaction)

        if reaction_template is None:
            self.ReactionTemplate = Classes.ReactionTemplate('none', '', [], [])
        else:
            self.ReactionTemplate = reaction_template
        self.Data = info

        self.Reactants = [Chem.MolToSmiles(x, canonical = True)
                                        for x in rdkit_reaction.GetReactants()]
        self.Products  = [Chem.MolToSmiles(x, canonical = True)
                                        for x in rdkit_reaction.GetProducts()]

