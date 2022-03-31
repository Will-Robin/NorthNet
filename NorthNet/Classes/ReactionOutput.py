class ReactionOutput:
    """
    Class to store outputs as reactions into reaction network
    """

    def __init__(self, reaction_output_string):
        """

        Parameters
        ----------
        reaction_input_string: str
            token for reaction input
            should follow the convention
            SMILES>>#0
        """
        assert isinstance(
            reaction_output_string, str
        ), """class ReactionOutput:
            reaction_output_string must be string like SMILES>>#0"""

        self.token = reaction_output_string

        token_sides = reaction_output_string.split(">>")
        self.OutputID = token_sides[1].split(".")
        self.OutputCompound = token_sides[0].split(".")

