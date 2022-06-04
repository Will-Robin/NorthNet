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

        assert isinstance(
            output, str
        ), """class NetworkOutput:
            id arg must be str like SMILES_#0."""

        self.token = output
        self.In = []
