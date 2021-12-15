import sys

class NetworkOutput:
    '''
    Class to store inputs into reaction network
    '''
    def __init__(self,output):
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

        if not isinstance(output, str):
            sys.exit('''class NetworkOutput:
                    id arg must be str like SMILES_#0.''')

        self.Mol = None

        # Note that the self.SMILES does not
        # actually hold a valid SMILES string,
        # it stores the input ID code
        self.SMILES = output

        self.In = []
        self.Out = None

        self.ReactiveSubstructures = None


