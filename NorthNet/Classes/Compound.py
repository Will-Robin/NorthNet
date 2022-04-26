from rdkit import Chem
import NorthNet.molecule_operations as mol_ops


class Compound:
    """
    Class to store compound information.
    """

    def __init__(self, SMILES):
        """
        Parameters
        ----------
        mol: str
            SMILES corresponding to compound.
        """

        assert isinstance(
            SMILES, str
        ), "class Compound: argument should be a SMILES string."

        self.Mol = Chem.MolFromSmiles(SMILES)

        canonical_smiles = mol_ops.canonicalise(SMILES)

        self.SMILES = canonical_smiles

        self.In = []
        self.Out = []

        self.ReactiveSubstructures = None
