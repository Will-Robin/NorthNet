import numpy as np
from NorthNet import data_processing as d_p
from NorthNet import info_params

class Dataset:
    '''
    Dataset object consisting of time series data for mutiple variables
    sharing the same time axis.
    '''
    def __init__(self, name, conditions, data, time_key = 'time/ s', get_flow_profile = True):

        '''
        Parameters
        ---------------
        name: str
            A name for the data series
        conditions: dict
            A dictionary of the conditions under which the data were collected
        data: dict
            Data organised into a dictionary with time and individual variables
            accessible by keys.

        Attributes
        ----------
        self.name: str
            A name for the data object
        self.conditions: dict
            The conditions under which the data were collected.
        self.time: numpy array
            An array of time values for the data.
        self.dependents: dict
            A dictionary of dependent varialbes accessible by their names as
            keys.
        '''

        self.name = name
        self.conditions = conditions

        if 'dilution_factor' in conditions:
            self.dilution = float(conditions['dilution_factor'][0])
            del self.conditions["dilution_factor"]
        if 'internal_ref_concentration/ M' in conditions:
            self.internal_ref_concentration = float(conditions['internal_ref_concentration/ M'][0])
            del self.conditions["internal_ref_concentration/ M"]

        self.time = data[time_key]
        self.independent = data[time_key]
        self.independent_name = time_key
        self.dependents = {k:data[k] for k in data if k != time_key}

        if get_flow_profile:
            self.get_input_profile()
        else:
            pass#print("flow profile not extracted (see Dataset object)")

    def get_input_profile(self):
        '''
        Adds flow profile attributes to the object using the conditions attribute.
        '''

        inds = np.where((np.array(self.conditions["flow_profile_time/ s"]) > np.amin(self.time))&(np.array(self.conditions["flow_profile_time/ s"]) < np.amax(self.time)))[0]

        self.input_flows_time = np.array(self.conditions["flow_profile_time/ s"])[inds]

        self.input_flows = {}

        for c in self.conditions:
            if "flow" in c and not "time" in c:
                profile = np.array(self.conditions[c])

                new_x, new_y = d_p.interpolate_traces(self.input_flows_time, profile[inds], length = len(self.time))

                self.input_flows[c] = new_y

        self.input_flows_time,_ = d_p.interpolate_traces(self.input_flows_time, self.input_flows_time, length = len(self.time))

        self.net_flow = np.zeros(len(self.input_flows[[*self.input_flows][0]]))
        for f in self.input_flows:
            self.net_flow = self.net_flow + self.input_flows[f]

    def get_mass_balance(self, inputs):

        '''
        flow_keys: list of tuples
            (c,d,e) c = concentration field, d = flow key, e = key to molecular_masses
        '''
        net_flow = self.net_flow/(60*60*10e6)
        reactor_vol = self.conditions["reactor_volume/ uL"][0]/10e6
        deltat = self.time[1]-self.time[0]

        mass_input = np.zeros(len(self.time))
        for c,d,e in inputs:
            mass_input = mass_input + info_params.molecular_masses[e]*self.conditions[c][0]*(self.input_flows[d]/(60*60*10e6))*deltat

        mass_in_reactor = np.zeros(len(self.time))

        mass_in_reactor[0] = reactor_vol*mass_input[0]/(deltat*net_flow[0])

        for x in range(0,len(mass_in_reactor)-1):
            mass_out = (mass_in_reactor[x]/reactor_vol)*net_flow[x]*deltat
            delta_mass = mass_input[x] - mass_out
            mass_in_reactor[x+1] = mass_in_reactor[x] + delta_mass

        self.mass_log = mass_in_reactor

        deps_mass = np.zeros(len(self.time))
        for d in self.dependents:
            deps_mass = deps_mass + self.dependents[d]*(self.conditions["reactor_volume/ uL"][0]/10e6)*info_params.molecular_masses[d[:-2]]

        self.measured_mass = deps_mass

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
        self.network = None
        self.path_scores = None

class DataReport:
    def __init__(self, file = ''):
        self.filename = 'not specified'
        self.experiment_code = 'not specified'
        self.conditions = {}
        self.analysis_details = {}
        self.series_values = np.array([])
        self.series_unit = 'not specified'
        self.data = {}

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
        file: path to file
        '''
        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']

        self.filename = str(file)

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

        e = [list(i) for i in zip(*dataset)]
        d_out = {}
        for s in e:
            d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

        self.series_unit = dataset[0][0]
        self.series_values = d_out[self.series_unit]

        del d_out[self.series_unit]

        self.data = d_out

        errors = self.import_file_section(file, "start_errors", "end_errors")

        e = [list(i) for i in zip(*dataset)]
        errors_out = {}
        for s in e:
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
        an_type = self.analysis_details['Chromatography_method'][0]

        if filename == '':
            filename = self.filename
        elif not filename.endswith('csv'):
            filename = filename + 'csv'
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
        entries = []
        repeat_entries = []
        for d in self.data:
            token = d.split(' ')[0]
            if token in entries:
                repeat_entries.append(token)
            entries.append(token)

        return list(set(repeat_entries))

    def remove_repeat_entries(self):
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
        Parameters
        ----------
        remove_list: list
            List of entries to remove from self.data
        '''
        for r in remove_list:
            del self.data[r]

    def remove_entries_below_threshold(self, threshold):
        '''
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
