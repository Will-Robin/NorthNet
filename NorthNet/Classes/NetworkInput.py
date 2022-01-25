from rdkit import Chem

class NetworkInput:
    '''
    Class to store inputs into reaction network
    '''
    def __init__(self,net_input):
        '''
        Designed to behave like a compound object
        self.SMILES does not actually hold a valid SMILES string,
        it stores the input ID code (e.g. SMILES_#0).

        Parameters
        ----------
        id: str
            token for reaction input
            should follow the convention
            SMILES_#0
        '''

        assert isinstance(net_input, str), \
            '''class NetworkInput:
            id arg must be str like SMILES_#0.'''

        self.Mol = Chem.MolFromSmiles(net_input.split('_')[0])

        # Note that the self.SMILES does not
        # actually hold a valid SMILES string,
        # it stores the input ID code
        self.SMILES = net_input

        self.In = None
        self.Out = []

        self.ReactiveSubstructures = None

