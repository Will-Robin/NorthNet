from rdkit.Chem import AllChem


class ReactionTemplate:
    """
    Class for reaction templates.
    """

    def __init__(self, name, reaction_SMARTS, reactant_substructs, product_substructs):
        """
        Parameters
        ----------
        name: str
            name for reaction.
        reaction_SMARTS: str
            reaction SMARTS string.

        reactant_substructs: list of strings
            reacting substructures

        product_substructs: list of strings
            substructures to which reactant substructures are converted.
        """

        assert isinstance(
            name, str
        ), """class ReactionTemplate:
            name arg should be a name for the ReactionTemplate as a string."""

        assert isinstance(
            reaction_SMARTS, str
        ), """class ReactionTemplate:
            the reaction_SMARTS arg should be a reaction SMARTS string."""

        if isinstance(reactant_substructs, list):
            check_reactant_substructs = [
                isinstance(r, str) for r in reactant_substructs
            ]
            assert all(
                check_reactant_substructs
            ), """class ReactionTemplate:
                reactant_substructs arg should be a list of SMARTS strings."""
        else:
            assert isinstance(
                reactant_substructs, list
            ), """class ReactionTemplate:
            the reactant_substructs arg should be a list of SMARTS strings."""

        if isinstance(product_substructs, list):
            check_product_substructs = [isinstance(s, str) for s in product_substructs]
            assert all(
                check_product_substructs
            ), """class ReactionTemplate:
                product_substructs arg should be a list of SMARTS strings."""
        else:
            assert isinstance(
                product_substructs, list
            ), """class ReactionTemplate:
                product_substructs arg should be a list of SMARTS strings."""

        self.Name = name
        if reaction_SMARTS != "":
            self.Reaction = AllChem.ReactionFromSmarts(reaction_SMARTS)
            self.ReactionSMARTS = AllChem.ReactionToSmarts(self.Reaction)
        else:
            self.Reaction = None
            self.ReactionSMARTS = ""

        self.ReactantSubstructures = reactant_substructs
        self.ProductSubstructures = product_substructs
