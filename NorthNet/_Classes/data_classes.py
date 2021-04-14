import numpy as np
from NorthNet import data_processing as d_p
from NorthNet.information import info_params

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
