import numpy as np


class DataSet:
    """
    Container for multiple data reports
    """

    def __init__(self, data_reports=[]):
        """
        data_reports: list of NorthNet DataReport objects
            data reports to create the data set
        """
        self.data_reports = dict()
        self.compounds = []

        if len(data_reports) == 0:
            pass
        else:
            for d in data_reports:
                self.add_data_report(d)

    def add_data_report(self, data_report):
        """
        Add data report to DataSet

        data_report: NorthNet DataReport
            DataReport to be added.
        """
        self.data_reports[len(self.data_reports) + 1] = data_report

        for d in data_report.data:
            if d not in self.compounds:
                self.compounds.append(d)

    def find_entry(self, entry_name):
        """
        Get the series_values and values of a given compound.

        entry_name: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        Returns: tuple of 1D numpy arrays
            (series_values, variable)
        """
        x_ax = np.array([])
        y_ax = np.array([])
        for d in self.data_reports:
            if entry_name in self.data_reports[d].data:
                y_ax = np.hstack((y_ax, self.data_reports[d].data[entry_name]))
                x_ax = np.hstack((x_ax, self.data_reports[d].series_values))

        return x_ax, y_ax

    def get_entry(self, entry_name):
        """
        Wrapper for find_entry(). Sorts the arrays so x_ax values are increasing
        with increasing index.

        entry_name: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        Returns: tuple of 1D numpy arrays
            (series_values, variable)
        """

        x_ax, y_ax = self.find_entry(entry_name)

        i = np.argsort(x_ax)
        x_ax = x_ax[i]
        y_ax = y_ax[i]

        return x_ax, y_ax

    def get_entry_indices(self, entry):
        """
        Get the indices which order the data so x_ax values are increasing
        with increasing index.

        entry: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        returns: i
            numpy array of int
        """
        x_ax, _ = self.find_entry(entry)
        i = np.argsort(x_ax)
        return i
