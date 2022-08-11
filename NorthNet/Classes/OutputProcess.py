class OutputProcess:
    """
    Class to store outputs as reactions into reaction network
    """

    def __init__(self, reaction_output_string):
        """

        Parameters
        ----------
        reaction_input_string: str
            A token for reaction output should follow the convention
            `SMILES>>#0`

        Attributes
        ----------
        token: string
            Token for the output.
        OutputID: str
            ID for the output.
        OutputCompound: str
            The compound associated with the output.
        In: list[str]
            List of processes connected to the output.
        """

        assert isinstance(
            reaction_output_string, str
        ), """class ReactionOutput:
            reaction_output_string must be string like SMILES>>#0"""

        self.token = reaction_output_string

        token_sides = reaction_output_string.split(">>")
        self.OutputID = token_sides[1]
        self.OutputCompound = token_sides[0]

        self.In = []
