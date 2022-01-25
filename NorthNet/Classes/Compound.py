from rdkit import Chem

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

        assert isinstance(SMILES, str), \
                'class Compound: argument should be a SMILES string.'

        self.Mol = Chem.MolFromSmiles(SMILES)
        if self.Mol is None:
            self.SMILES = SMILES
        else:
            self.SMILES = Chem.MolToSmiles(self.Mol)

        self.In = []
        self.Out = []

        self.ReactiveSubstructures = None

