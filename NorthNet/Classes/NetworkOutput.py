
class NetworkOutput:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, output):
        """

        Parameters
        ----------
        id: str
            token for reaction input
            should follow the convention
            SMILES_#0
        """

        assert isinstance(
            output, str
        ), """class NetworkOutput:
            id arg must be str like SMILES_#0."""

        self.token = output
        self.In = []

