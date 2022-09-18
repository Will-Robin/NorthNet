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

        self.verify_input()

        self.token = net_input
        self.Out = []

    def verify_input(self, net_input):
        """
        Check that the input arguments are of the correct type.

        Parameters
        ----------
        net_input: str
        """

        assert isinstance(
            net_input, str
        ), """class NetworkInput:
            id arg must be str like SMILES_#0."""
