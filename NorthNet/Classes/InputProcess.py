class InputProcess:
    """
    Class to store inputs as reactions into reaction network
    """

    def __init__(self, reaction_input_string):
        """

        Parameters
        ----------
        reaction_input_string: str
            A token for reaction input should follow the convention
            `SMILES_#0>>SMILES`

        Attributes
        ----------
        InputID: str
            ID for the input.
        InputCompound: str
            Compound parameters associated with the input.
        Out: list[str]
            List of compounds arising from the input.
        """

        self.verify_input(reaction_input_string)

        self.token = reaction_input_string

        token_sides = reaction_input_string.split(">>")
        self.InputID = token_sides[0]
        self.InputCompound = token_sides[1].split(".")

        self.Out = []

    def verify_input(self, reaction_input_string):
        """
        Check that the input arguments are of the correct type.

        Parameters
        ----------
        reaction_input_string: str
        """

        assert isinstance(
            reaction_input_string, str
        ), """class ReactionInput:
            reaction_input_string arg must be str like SMILES_#0>>SMILES."""
