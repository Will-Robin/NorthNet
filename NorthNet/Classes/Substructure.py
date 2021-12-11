import sys
from rdkit import Chem

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

        if not isinstance(SMARTS, str):
            sys.exit('class Substructure: argument should be a SMARTS string.')

        self.Mol = Chem.MolFromSmarts(SMARTS)

        if Chem.MolToSmarts(self.Mol) == None:
            self.SMARTS = SMARTS
        else:
            self.SMARTS = Chem.MolToSmarts(self.Mol)

        self.In = []
        self.Out = []

