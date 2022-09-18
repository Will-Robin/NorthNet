class NetworkOutput:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, output):
        """

        Parameters
        ----------
        output: str
            A token for reaction input should follow the convention `SMILES_#0`

        Attributes
        ----------
        token: str

        In: list[str]
            List of Reaction Output processes connected to the output.
        """

        self.verify_input(output)

        self.token = output
        self.In = []

    def verify_input(self, output):
        """
        Check that the input arguments are of the correct type.

        Parameters
        ----------
        output: str
        """

        assert isinstance(
            output, str
        ), """class NetworkOutput:
            id arg must be str like SMILES_#0."""
