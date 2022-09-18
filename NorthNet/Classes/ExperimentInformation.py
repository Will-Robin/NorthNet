class ExperimentInformation:
    def __init__(self, name, path, parameters, modulation):
        """

        Parameters
        ----------
        name: str
        path: str
        parameters: dict
        modulation: str

        Attributes
        ----------
        name: str
        path: str
        parameters: dict
        modulation: str
        """
        self.name = name
        self.path = path
        self.parameters = parameters
        self.modulation = modulation
