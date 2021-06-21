import numpy as np
from pathlib import Path

class Experiment_Information:
    def __init__(self, name, path, parameters, modulation):
        '''
        name: str
        path: str
        parameters: dict
        '''
        self.name = name
        self.path = path
        self.parameters = parameters
        self.modulation = modulation

class DataReport:
    def __init__(self, file = ''):
        '''
        file: pathlib Path or str
            Path to file
        '''
        self.filename = 'not specified'
        self.experiment_code = 'not specified'
        self.conditions = {}
        self.analysis_details = {}
        self.series_values = np.array([])
        self.series_unit = 'not specified'
        self.data = {}
        self.errors = {}

        if file == '':
            pass
        else:
            self.read_from_file(file)

    def import_file_section(self, file, start_token, end_token):
        '''
        Parameters
        ----------
        file: path to file
        start_token: str
            String in line to start reading file from.
        end_token:
            String in line to end reading file from.
        '''

        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        readstate = False
        c_set = []
        with open(file, 'r', encoding = 'latin-1') as f:
            for c,line in enumerate(f):
                if start_token in line:
                    readstate = True
                    line = next(f)
                if end_token in line:
                    readstate = False
                if readstate:
                    newline = spl_lin(line)
                    c_set.append(newline)

        return c_set

    def read_from_file(self, file):
        '''
        Parameters
        ----------
        file: pathlib Path or str
        '''
        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']

        if type(file) == str:
            file = Path(file)

        self.filename = file.name

        with open(file, 'r', encoding = 'latin-1') as f:
            for line in f:
                ins = spl_lin(line)
                if 'Dataset' in line:
                    self.experiment_code = ins[1]

        condset = self.import_file_section(file, "start_conditions",
                                           "end_conditions")

        c_out = {}
        for c in condset:
            entry = [float(x) for x in c[1:]]
            if len(entry) == 1:
                self.conditions[c[0]] = entry[0]
            else:
                self.conditions[c[0]] = np.array(entry)

        dataset = self.import_file_section(file, "start_data", "end_data")

        transposed_datalines = [list(i) for i in zip(*dataset)]
        d_out = {}
        for s in transposed_datalines:
            d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

        self.series_unit = dataset[0][0]
        self.series_values = d_out[self.series_unit]

        del d_out[self.series_unit]

        self.data = d_out

        errors = self.import_file_section(file, "start_errors", "end_errors")

        if len(errors) == 0:
            self.errors = {d:np.zeros(len(self.series_values)) for d in self.data}
        else:
            transposed_error_lines = [list(i) for i in zip(*errors)]
            errors_out = {}
            for s in transposed_error_lines:
                errors_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])
            del errors_out[self.series_unit]
            self.errors = errors_out

        analysis = self.import_file_section(file, "start_analysis_details",
                                                         "end_analysis_details")
        for a in analysis:
            self.analysis_details[a[0]] = [x for x in a[1:]]

    def write_conditions_header(self, outfile):
        '''
        Parameters
        ----------
        outfile: Python file object
        '''
        # writing experiment conditions to file
        outfile.write("Dataset,{}\n".format(self.experiment_code))
        outfile.write("start_conditions\n")
        for c in self.conditions:
            outfile.write("{},".format(c))
            if type(self.conditions[c]) == float:
                outfile.write("{},".format(self.conditions[c]))
            else:
                [outfile.write("{},".format(x)) for x in self.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")
        # writing analysis details
        outfile.write("start_analysis_details\n")
        for ad in self.analysis_details:
            outfile.write('{},'.format(ad))
            if type(self.analysis_details[ad]) == str:
                outfile.write('{},'.format(self.analysis_details[ad]))
            else:
                [outfile.write('{},'.format(x)) for x in self.analysis_details[ad]]
            outfile.write('\n')
        outfile.write("end_analysis_details\n")

    def write_to_file(self, filename = '', path = None):
        '''
        Parameters
        ----------
        filename: str
            name for file
        path: pathlib Path object
            Path to folder for file storage.
        '''

        import numpy as np
        if 'Chromatography_method' in self.analysis_details:
            an_type = self.analysis_details['Chromatography_method'][0]
        else:
            an_type = 'not specified'

        if filename == '':
            filename = self.filename
        elif not filename.endswith('.csv'):
            filename = filename + '.csv'
        if path == None:
            fname = filename
        else:
            fname = path/filename

        with open(fname, 'w') as outfile:
            # writing experiment conditions to file
            self.write_conditions_header(outfile)
            # writing data
            sorted_keys = sorted([*self.data], key = lambda x:x.count('C'))

            outfile.write("start_data\n")

            p_header = [self.series_unit]
            out = np.array([self.series_values])

            for s in sorted_keys:
                p_header.append(s)
                out = np.vstack((out,self.data[s]))

            out = out.T
            [outfile.write("{},".format(x)) for x in p_header]

            outfile.write("\n")

            for x in range(0,len(out)):
                for y in range(0,len(out[x])):
                    outfile.write("{},".format(out[x,y]))
                outfile.write("\n")

            outfile.write("end_data\n")

            outfile.write('start_errors\n')
            if len(self.errors)> 0:
                err_out = np.array([self.series_values])

                for s in sorted_keys:
                    err_out = np.vstack((err_out,self.errors[s]))

                err_out = err_out.T
                [outfile.write("{},".format(x)) for x in p_header]
                outfile.write("\n")
                for x in range(0,len(err_out)):
                    for y in range(0,len(err_out[x])):
                        outfile.write("{},".format(err_out[x,y]))
                    outfile.write("\n")

            outfile.write('end_errors\n')

    def find_repeat_data_entries(self):
        '''
        Find compound entries which are repeated in self.data
        '''
        entries = []
        repeat_entries = []
        for d in self.data:
            token = d.split(' ')[0]
            if token in entries:
                repeat_entries.append(token)
            entries.append(token)

        return list(set(repeat_entries))

    def remove_repeat_entries(self):
        '''
        Remove compound entries which are repeated in self.data
        '''
        import numpy as np
        # deleting duplicate entries: taking the entry with the higher signal using the
        # signal sum as a discriminant.
        repeat_entries = self.find_repeat_data_entries()

        for r in repeat_entries:
            compare_keys = []
            for d in self.data:
                if r in d:
                    compare_keys.append(d)

            checkline = np.zeros(len(compare_keys))
            for c,comp in enumerate(compare_keys):
                checkline[c] = np.sum(self.data[comp])

            i_min = np.argmin(checkline)

            del self.data[compare_keys[i_min]]

    def remove_specific_entries(self,remove_list):
        '''
        remove entries in remove list from self.data

        Parameters
        ----------
        remove_list: list
            List of entries to remove from self.data
        '''
        for r in remove_list:
            del self.data[r]

    def remove_entries_below_threshold(self, threshold):
        '''
        remove entries whose maximum value does not exceed threshold

        Parameters
        ----------
        threshold: float
            threshold below which entries will be removed.
        '''
        import numpy as np
        # remove entries whose concentrations/integrals do not cross a defined boundary
        del_list = []
        for d in self.data:
            if np.amax(self.data[d]) < threshold:
                del_list.append(d)

        self.remove_specific_entries(del_list)

class DataSet:
    '''
    Container for multiple data reports
    '''
    def __init__(self, data_reports = []):
        '''
        data_reports: list of NorthNet DataReport objects
            data reports to create the data set
        '''
        self.data_reports = {}
        self.compounds = []

        if len(data_reports) == 0:
            pass
        else:
            for d in data_reports:
                self.add_data_report(d)


    def add_data_report(self, data_report):
        '''
        Add data report to DataSet

        data_report: NorthNet DataReport
            DataReport to be added.
        '''
        self.data_reports[len(self.data_reports)+1] = data_report

        for d in data_report.data:
            if d not in self.compounds:
                self.compounds.append(d)

    def find_entry(self,entry_name):
        '''
        Get the series_values and values of a given compound.

        entry_name: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        Returns: tuple of 1D numpy arrays
            (series_values, variable)
        '''
        x_ax = np.array([])
        y_ax = np.array([])
        for d in self.data_reports:
            if entry_name in self.data_reports[d].data:
                y_ax = np.hstack((y_ax, self.data_reports[d].data[entry_name]))
                x_ax = np.hstack((x_ax, self.data_reports[d].series_values))

        return x_ax, y_ax

    def get_entry(self, entry_name):
        '''
        Wrapper for find_entry(). Sorts the arrays so x_ax values are increasing
        with increasing index.

        entry_name: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        Returns: tuple of 1D numpy arrays
            (series_values, variable)
        '''

        x_ax, y_ax = self.find_entry(entry_name)

        i = np.argsort(x_ax)
        x_ax = x_ax[i]
        y_ax = y_ax[i]

        return x_ax, y_ax

    def get_entry_indices(self,entry):
        '''
        Get the indices which order the data so x_ax values are increasing
        with increasing index.

        entry: str
            Key to the compound in DataReport.data, e.g. [C=O]/ M

        returns: i
            numpy array of int
        '''
        x_ax, y_ax = self.find_entry(entry)
        i = np.argsort(x_ax)
        return i
