from rdkit import Chem
import NorthNet.molecule_operations as mol_ops


class Compound:
    """
    A class to store compound information.
    """

    def __init__(self, SMILES):
        """
        Parameters
        ----------
        SMILES: str
            SMILES corresponding to compound.

        Attributes
        ----------
        Mol: rdkit.Chem.rdchem.Mol
            RDKit representation of the compound to allow for chemoinformatics
            operations.
        In: list[str]
            List of reaction tokens for which the compound is a product.
        Out: list[str]
            List of reaction tokens for which the compound is a reactant.
        ReactiveSubstructures: list[str]
            List of keys to contextual reactive substructures which either
            engage in, or are products of, reactions.
        """

        assert isinstance(
            SMILES, str
        ), "class Compound: argument should be a SMILES string."

        self.Mol = Chem.MolFromSmiles(SMILES)

        canonical_smiles = mol_ops.canonicalise(SMILES)

        self.SMILES = canonical_smiles

        self.In = []
        self.Out = []

        self.ReactiveSubstructures = []
