from rdkit import Chem


class Substructure:
    """
    Class to store substructures.
    """

    def __init__(self, SMARTS):
        """
        Parameters
        ----------
        SMARTS: str
            SMARTS corresponding to substructure.
        """

        assert isinstance(
            SMARTS, str
        ), """class Substructure:
            argument should be a SMARTS string."""

        self.Mol = Chem.MolFromSmarts(SMARTS)

        if Chem.MolToSmarts(self.Mol) is None:
            self.SMARTS = SMARTS
        else:
            self.SMARTS = Chem.MolToSmarts(self.Mol)

        self.In = []
        self.Out = []
