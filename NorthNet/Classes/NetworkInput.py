class NetworkInput:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, net_input):
        """

        Parameters
        ----------
        net_input: str
            A token for reaction input. It should follow the convention
            `SMILES_#0`.

        Attributes
        ----------
        token: str
        Out: list[str]
            List of reaction input tokens which this NetworkInput feeds into.
        """

        assert isinstance(
            net_input, str
        ), """class NetworkInput:
            id arg must be str like SMILES_#0."""

        self.token = net_input
        self.Out = []
