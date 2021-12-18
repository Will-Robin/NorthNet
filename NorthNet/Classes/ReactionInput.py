
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

        assert isinstance(id, str), \
            '''class ReactionInput:
            reaction_input_string arg must be str like SMILES_#0>>SMILES.'''

        self.Reaction = None
        self.ReactionSMILES = reaction_input_string
        self.Reactants, self.Products = self.reactants_products_from_string(
                                                        reaction_input_string)

        self.CompoundInput = self.Products[0]
        self.InputID = self.Reactants[0]

        self.ReactionTemplate = None
        self.Data = {}

    def reactants_products_from_string(self, reaction_smiles):
        '''
        C
        '''
        split_rxn_smiles = reaction_smiles.split('>>')
        reactants = [x for x in split_rxn_smiles[0].split('.') if x != '']
        products = [x for x in split_rxn_smiles[1].split('.') if x != '']
        return reactants, products

