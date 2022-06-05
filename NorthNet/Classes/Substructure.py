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

        Attributes
        ----------
        Mol: rdkit.Chem.rdchem.Mol
        SMARTS: str
        MatchingCompounds: list[str]
            Contextual list of compound tokens which contain this substructure.
        ReactionParticipations: list[str]
            List of reaction tokens for which the substructure is a reactant.
        ProducingReactions: list[str]
            List of reaction tokens of which the substructure is a product.
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

        self.MatchingCompounds = []
        self.ReactionParticipations = []
        self.ProducingReactions = []
