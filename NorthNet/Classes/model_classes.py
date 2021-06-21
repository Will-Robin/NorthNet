class ModelWriter:
    def __init__(self, network = None, experiment = None,
                       input_token = '_#0',
                       output_token = 'Sample',
                       flowrate_time_conversion = 3600,
                       time_limit = False,
                       lead_time = 1000):
        '''
        network: NorthNet Network
        experiment: NorthNet DataReport,
        input_token: str
        output_token: str
        flowrate_time_conversion: float
            Conversion for the time component:
            CAUTION: this class currently expects flow rates to be given
            in units of uL/h and reactor volumes to be given in uL,
            so conversion errors may result if the input DataReport's
            attributes fall out of this pattern.
        time_limt: bool or float
            How far in time the flow profile will be considered
            in generating the model.
        '''

        self.network = network
        self.input_token = input_token
        self.output_token = output_token
        self.flowrate_time_conversion = flowrate_time_conversion
        self.time = np.inf
        self.time_limit = time_limit
        self.lead_time = lead_time
        self.observed_compounds = []
        self.reactor_volume = 1.0
        self.flow_profile_time = []
        self.flow_profiles = {}
        self.sigma_flow = []
        self.name = ''
        self.species = {}
        self.rate_constants = {}
        self.inputs = {}
        self.inflows = {}
        self.outflows = {}
        self.time_offset = 0.0

        if network == None:
            pass
        else:
            self.get_network_tokens()

        if experiment == None:
            pass
        else:
            self.load_experiment_details(experiment)

    def get_network_tokens(self):
        '''

        '''
        network = self.network
        compounds = [network.NetworkCompounds[x] for x in network.NetworkCompounds]
        network_inputs = [*network.NetworkInputs]
        inflows = [r for r in network.NetworkReactions if self.input_token in r]
        outflows = [r for r in network.NetworkReactions if self.output_token in r]
        reactions = [r for r in network.NetworkReactions]

        for i in inflows:
            reactions.remove(i)
        for o in outflows:
            reactions.remove(o)

        species = {s.SMILES:"S[{}]".format(c) for c,s in enumerate(compounds) if s != ''}
        rate_consts = {k:"k[{}]".format(c) for c,k in enumerate(reactions)}
        inputs = {i:0.0 for c,i in enumerate(network_inputs)}
        flow_ins = {i:'I[{}]'.format(c) for c,i in enumerate(inflows)}
        flow_outs = {o:'sigma_flow' for o in outflows}

        self.name = network.Name
        self.species = species
        self.rate_constants = rate_consts
        self.inputs = inputs
        self.inflows = flow_ins
        self.outflows = flow_outs

    def load_experiment_details(self, experiment):
        if experiment.series_unit == 'time/ s':
            self.time = experiment.series_values.copy()

        for d in experiment.data:
            self.observed_compounds.append(d.split('/')[0].split(' ')[0])

        for c in experiment.conditions:
            if 'reactor_volume' in c:
                self.reactor_volume = experiment.conditions[c]
            elif ' M' in c:
                standardised_key = c.split('_')[0].split('/')[0]
                clef = info_params.canonical_SMILES[standardised_key]
                for f in self.inputs:
                    stand_flow_key = f.split('_')[0]
                    if clef == stand_flow_key:
                        self.inputs[f] = experiment.conditions[c]
            elif 'time' in c and 'flow' in c:
                self.flow_profile_time = experiment.conditions[c].copy()

            elif 'flow' in c and not 'time' in c:
                standardised_key = c.split(' ')
                clef = info_params.canonical_SMILES[standardised_key[0].split('_')[0]]
                self.flow_profiles[clef] = experiment.conditions[c].copy()

        if self.time_limit:
            t_lim_max = min(np.amax(self.time), self.time_limit)
        else:
            t_lim_max = np.amax(self.time)

        t_lim_min = self.time[0] - self.lead_time

        idx = np.where((self.flow_profile_time > t_lim_min)&
                       (self.flow_profile_time < t_lim_max))[0]

        self.flow_profile_time = self.flow_profile_time[idx]
        self.time_offset = self.flow_profile_time[0]
        self.flow_profile_time -= self.time_offset
        self.time -= self.time_offset

        self.sigma_flow = np.zeros(len(self.flow_profile_time))
        for f in self.flow_profiles:
            self.flow_profiles[f] = self.flow_profiles[f][idx]/self.reactor_volume
            self.flow_profiles[f] /= self.flowrate_time_conversion
            self.sigma_flow += self.flow_profiles[f]

    def write_equation_text(self):

        '''
        Parameters
        ----------

        Returns
        -------
        eq_text: str
            Rate equations in text form.
        '''

        network = self.network
        compounds = [x for x in network.NetworkCompounds]
        reactions = [*network.NetworkReactions]

        eq_text = ""

        for count,c in enumerate(compounds):
            eq_text += "P[{}] = ".format(count)
            for i in network.NetworkCompounds[c].In:
                if '_#0' in i:
                    in_compound = network.NetworkReactions[i].InputID
                    ki = '+{}*{}'.format(self.inflows[i], self.inputs[in_compound])
                    eq_text += ki
                else:
                    reactants = network.NetworkReactions[i].Reactants
                    # remove water from reactants
                    reactants = [x for x in reactants if x != 'O']

                    ki = "+{}*".format(self.rate_constants[i])

                    if len(reactants) == 0:
                        specs = ''#inflows[i]
                    else:
                        specs = "*".join([self.species[x] for x in reactants])

                    eq_text += "{}{}".format(ki,specs)

            for o in network.NetworkCompounds[c].Out:
                if 'Sample' in o:
                    out_compound = network.NetworkReactions[o].CompoundOutput
                    ki = '-{}*{}'.format(self.outflows[o], self.species[out_compound])
                    eq_text += ki
                else:
                    ki = "-{}*".format(self.rate_constants[o])
                    specs = "*".join([self.species[x] for x in network.NetworkReactions[o].Reactants])
                    eq_text += "{}{}".format(ki,specs)

            eq_text += "\n"

        return eq_text

    def write_flow_profile_text(self):
        collection_array = np.zeros((len(self.flow_profiles)+2,
                                     len(self.flow_profile_time)))

        collection_array[0] = self.flow_profile_time
        for c,f in enumerate(self.inputs,1):
            collection_array[c] = self.flow_profiles[f.split('_')[0]]

        collection_array[-1] = self.sigma_flow

        text = 'F = np.array('
        text += np.array2string(collection_array,
                                     formatter={'float_kind':lambda x: "%.9f" % x},
                                     separator=',',threshold=np.inf)
        text += ')'

        return text

    def write_model_equation_text(self):
        '''
        Returns
        -------
        eq_text: str
            Rate equations in text form.
        '''
        network = self.network

        compounds = [x for x in network.NetworkCompounds]
        reactions = [*network.NetworkReactions]


        eq_text = ""

        for count,c in enumerate(compounds):
            eq_text += "P[{}] = ".format(count)
            for i in network.NetworkCompounds[c].In:
                if '_#0' in i:
                    in_compound = network.NetworkReactions[i].InputID
                    ki = '+{}*{}'.format(self.inflows[i], self.inputs[in_compound])
                    eq_text += ki
                else:
                    reactants = network.NetworkReactions[i].Reactants
                    # remove water from reactants
                    reactants = [x for x in reactants if x != 'O']

                    ki = "+{}*".format(self.rate_constants[i])

                    if len(reactants) == 0:
                        specs = ''#inflows[i]
                    else:
                        specs = "*".join([self.species[x] for x in reactants])

                    eq_text += "{}{}".format(ki,specs)

            for o in network.NetworkCompounds[c].Out:
                if 'Sample' in o:
                    out_compound = network.NetworkReactions[o].CompoundOutput
                    ki = '-{}*{}'.format(self.outflows[o], self.species[out_compound])
                    eq_text += ki
                else:
                    ki = "-{}*".format(self.rate_constants[o])
                    specs = "*".join([self.species[x] for x in network.NetworkReactions[o].Reactants])
                    eq_text += "{}{}".format(ki,specs)

            eq_text += "\n"

        return eq_text

    def write_to_module_text_A(self, numba_decoration = False):
        get_index = lambda x: int(x[x.find("[")+1:x.find("]")])

        flow_profile_text = self.write_flow_profile_text()
        model_text = self.write_model_equation_text().split('\n')

        lines = ["import numpy as np\n"]
        if numba_decoration:
            lines.append("import numba\n\n")
            lines.append("@numba.jit(numba.float64[:](numba.float64,numba.float64[:],numba.float64[:]),\n"
                       "\tlocals={'P': numba.float64[:],'F': numba.float64[:,:],'I':numba.float64[:]},nopython=True)\n")
        lines.append("def model_function(time, S, k):\n\n")
        lines.append("\tP = np.zeros(len(S))\n\n")
        lines.append("\t")
        lines.append(flow_profile_text)
        lines.append("\n")
        lines.append("\n")
        lines.append("\tidx = np.abs(F[0] - time).argmin()\n")
        lines.append("\n")
        lines.append("\tI = F[1:-1,idx]\n")
        lines.append("\n")
        lines.append('\tsigma_flow = F[-1,idx]\n')
        lines.append("\n")
        lines.append("\n")

        for m in model_text:
            lines.append("\t{}\n".format(m))

        lines.append("\treturn P\n")
        lines.append("\n")
        lines.append("def wrapper_function(time, S, k):\n")
        lines.append("\treturn model_function(time, S, k)\n")
        lines.append("\n")

        lines.append("\n")
        lines.append("species = {")
        for k in self.species:
            idx = get_index(self.species[k])
            lines.append("'{}':{},".format(k,idx))
        lines.append("}\n")

        lines.append("\n")
        lines.append("reactions = {")
        for k in self.rate_constants:
            idx = get_index(self.rate_constants[k])
            lines.append("'{}':{},".format(k,idx))
        lines.append("}\n")

        lines.append("\n")
        lines.append("inputs = {")
        for k in self.inputs:
            idx = self.inputs[k]
            lines.append("'{}':{},".format(k,idx))

        lines.append("}\n")
        lines.append("\n")

        lines.append("k = np.zeros(max(reactions.values())+1) # rate constants\n")
        lines.append("\n")
        lines.append("S = np.zeros(len(species)) # initial concentrations\n")
        lines.append("\n")

        lines.append("C = np.zeros(len(inputs)) # input concentrations\n")
        lines.append("\n")
        lines.append("time_offset = {}\n".format(self.time_offset))
        lines.append("lead_in_time = {}\n".format(self.lead_time))

        text = ''
        for l in lines:
            text += l

        return text