class NetworkInput:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, net_input):
        """

        Parameters
        ----------
        id: str
            token for reaction input
            should follow the convention
            SMILES_#0
        """

        assert isinstance(
            net_input, str
        ), """class NetworkInput:
            id arg must be str like SMILES_#0."""

        self.token = net_input
        self.Out = []
