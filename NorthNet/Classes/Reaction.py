from NorthNet import Classes
from NorthNet.text_parsing import conversions as conv
from NorthNet import reaction_operations as reac_ops


class Reaction:
    """
    A class representing chemical reactions
    """

    def __init__(self, reaction_smiles, reaction_template=None, info=dict()):
        """
        Parameters
        ----------
        reaction_smiles: str
            Reaction object for the reaction. Please provide a valid reaction
            SMILES string with canonicalised SMILES reactants and products.
        reaction_template: NorthNet Reaction_Template object
            Reaction template which created the reaction.
        info: dict
            Dictionary of information (e.g. database entries)

        Attributes
        ----------
        self.ReactionSMILES: str
            Reaction SMILES string
        self.Reaction_Template: NorthNet ReactionTemplate or None
            Reaction template for the reaction. Defaults to None
        self.info: dict
            Dictionary of information (e.g. database entries)
        """

        assert isinstance(
            reaction_smiles, str
        ), """class Reaction:
            reaction_smiles arg should be an string."""

        if reaction_template is not None:
            assert isinstance(
                reaction_template, Classes.ReactionTemplate
            ), """class Reaction:
            reaction_template should be None or a NortNet
            ReactionTemplate object."""

        assert isinstance(info, dict), """class Reaction: info arg should be a dict."""

        canonical_rxn_smiles = reac_ops.canonicalise_reaction_smiles(reaction_smiles)

        self.ReactionSMILES = canonical_rxn_smiles

        if reaction_template is None:
            self.ReactionTemplate = Classes.ReactionTemplate("none", "", [], [])
        else:
            self.ReactionTemplate = reaction_template

        self.Data = info

        self.Reactants, self.Products = conv.reaction_smiles_split(canonical_rxn_smiles)
