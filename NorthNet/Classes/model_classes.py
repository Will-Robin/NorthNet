import numpy as np

class ModelWriter:
    def __init__(self, network = None, experiment = None,
                       input_token = '_#0',
                       output_token = 'Sample',
                       flowrate_time_conversion = 3600,
                       time_limit = False,
                       lead_time = 1000):
        '''

        A class designed to generate modelling apparatus by combining a Network
        structure and experimental conditions.

        All of the input variables will be compiled into the model (including
        flow profile information). This may make a large object, but it should
        lighten the load in performing calculations.

        network: NorthNet Network
            
        experiment: NorthNet DataReport

        input_token: str

        output_token: str

        flowrate_time_conversion: float
            Conversion for the time component:

            CAUTION: this class currently expects flow rates to be given
            in units of uL/h and reactor volumes to be given in uL,
            so conversion errors may result if the input DataReport's
            attributes fall out of this pattern.

        time_limit: bool or float
            Time cutoff for the model calculation. The values of input flow
            profiles will be included below the time limit.

        lead_time: 1000
            Governs the time from which the model will begin calculation
            before the first data time point. (start time = t0 - lead_time)
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
            self.name = network.Name
            self.get_network_tokens()

        if experiment == None:
            pass
        else:
            self.load_experiment_details(experiment)

    def get_network_tokens(self):
        '''
        Get dictionaries of tokens for the compounds, reactions, inputs,
        inflows, outflows
        '''
        network = self.network
        compounds = [*network.NetworkCompounds]
        reactions = [*network.NetworkReactions]
        network_inputs = [*network.NetworkInputs]
        inflows = [r for r in reactions if self.input_token in r]
        outflows = [r for r in reactions if self.output_token in r]

        for i in inflows:
            reactions.remove(i)
        for o in outflows:
            reactions.remove(o)

        species = {s:f"S[{c}]" for c,s in enumerate(compounds) if s != ''}
        rate_consts = {k:f"k[{c}]" for c,k in enumerate(reactions)}
        inputs = {i:0.0 for c,i in enumerate(network_inputs)}
        flow_ins = {i:f'I[{c}]' for c,i in enumerate(inflows)}
        flow_outs = {o:'sigma_flow' for o in outflows}

        self.species = species
        self.rate_constants = rate_consts
        self.inputs = inputs
        self.inflows = flow_ins
        self.outflows = flow_outs

    def load_experiment_details(self, experiment):
        '''
        Load experiment details into the class attributes
        '''
        if experiment.series_unit == 'time/ s':
            self.time = experiment.series_values.copy()

        for d in experiment.data:
            self.observed_compounds.append(d.split('/')[0].split(' ')[0])

        for c in experiment.conditions:
            if 'reactor_volume' in c:
                self.reactor_volume = experiment.conditions[c]
            elif '/ M' in c:
                smiles = c.split('/')[0]
                clef = smiles
                for f in self.inputs:
                    stand_flow_key = f.split('_')[0]
                    if smiles == stand_flow_key:
                        self.inputs[f] = experiment.conditions[c]
            elif 'time' in c and 'flow' in c:
                self.flow_profile_time = experiment.conditions[c].copy()

            elif 'flow' in c and not 'time' in c:
                smiles = c.split('_')[0]
                self.flow_profiles[smiles] = experiment.conditions[c].copy()

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

    def write_flow_profile_text(self, suffix = ""):
        '''
        Write flow profiles as a text numpy array.
        '''

        collection_array = np.zeros((len(self.flow_profiles)+2,
                                     len(self.flow_profile_time)))

        collection_array[0] = self.flow_profile_time
        for c,f in enumerate(self.inputs,1):
            collection_array[c] = self.flow_profiles[f.split('_')[0]]

        collection_array[-1] = self.sigma_flow

        text = 'F = np.array(\n'
        text += suffix
        text += np.array2string(collection_array,
                             formatter={'float_kind':lambda x: "%.9f" % x},
                             separator=',',
                             threshold=np.inf
                             ).replace("\n",f"\n{suffix}")
        text += ')'

        return text

    def write_model_equation_text(self):
        '''
        Writes models as equations with variables which refer to indices of
        arrays:

        P: 1D array of len(self.network.NetworkCompounds)
            Stores the product state if the system following calculation
        S: 1D array of len(self.network.NetworkCompounds)
            Stores the initial state of the system following calculation
            (see self.get_network_tokens())
        k: 1D array of len(self.network.NetworkReactions)
            rate constants arranged in standardised order.
            (see self.get_network_tokens())
        I: 1D array of len(self.network.NetworkInputs)
            Inputs into the system (see self.get_network_tokens())

        Includes output flow terms for all compounds.

        Returns
        -------
        eq_lines: list
            List of rate equations in text form.
        '''
        network = self.network

        compounds = [x for x in network.NetworkCompounds]
        reactions = [*network.NetworkReactions]

        eq_lines = []

        for count,c in enumerate(compounds):
            line_text = f"P[{count}] = "
            for i in network.NetworkCompounds[c].In:
                if '_#0' in i:
                    in_compound = network.NetworkReactions[i].InputID
                    ki = f'+{self.inflows[i]}*{self.inputs[in_compound]}'
                    line_text += ki
                else:
                    reactants = network.NetworkReactions[i].Reactants
                    # remove water from reactants
                    reactants = [x for x in reactants if x != 'O']

                    ki = f"+{self.rate_constants[i]}*"

                    if len(reactants) == 0:
                        specs = ''#inflows[i]
                    else:
                        specs = "*".join([self.species[x] for x in reactants])

                    line_text += f"{ki}{specs}"

            for o in network.NetworkCompounds[c].Out:
                if 'Sample' in o:
                    out_compound = network.NetworkReactions[o].CompoundOutput
                    ki = f'-{self.outflows[o]}*{self.species[out_compound]}'
                    line_text += ki
                else:
                    ki = f"-{self.rate_constants[o]}*"
                    reactants = [self.species[x]
                                for x in network.NetworkReactions[o].Reactants]
                    specs = "*".join(reactants)
                    line_text += f"{ki}{specs}"

            eq_lines.append(line_text)

        return eq_lines

    def write_variables_text(self):
        '''
        Write model variables as strings stored in a list
        '''
        get_index = lambda x: int(x[x.find("[")+1:x.find("]")])
        lines = []

        lines.append("species = {")
        for k in self.species:
            idx = get_index(self.species[k])
            lines.append(f"'{k}':{idx},")
        lines.append("}")

        lines.append("")
        lines.append("reactions = {")
        for k in self.rate_constants:
            idx = get_index(self.rate_constants[k])
            lines.append(f"'{k}':{idx},")
        lines.append("}")

        lines.append("")
        lines.append("inputs = {")
        for k in self.inputs:
            idx = self.inputs[k]
            lines.append(f"'{k}':{idx},")
        lines.append("}")
        lines.append("")

        lines.append("k = np.zeros(max(reactions.values())+1) # rate constants")
        lines.append("")
        lines.append("S = np.zeros(len(species)) # initial concentrations")
        lines.append("")

        lines.append("C = np.zeros(len(inputs)) # input concentrations")
        lines.append("")
        lines.append(f"time_offset = {self.time_offset}")
        lines.append(f"lead_in_time = {self.lead_time}")

        return lines

    def write_to_module_text(self, numba_decoration = False):
        '''
        Writing a Python script based on the object attributes
        This method creates an importable set of functions and varaibles
        which can be using in model calculations.

        Parameters
        ----------
        numba_decoration: bool
            Whether to include numbda compilation in the function.

        Returns
        -------
        text: str
            The module text.
        '''

        flow_profile_text = self.write_flow_profile_text(suffix = "\t\t")
        model_text = self.write_model_equation_text()

        lines = ["import numpy as np"]
        if numba_decoration:
            lines.append("import numba\n")
            lines.append("")
            lines.append("@numba.jit(numba.float64[:](numba.float64,numba.float64[:],numba.float64[:]),\n"
                       "\tlocals={'P': numba.float64[:],'F': numba.float64[:,:],'I':numba.float64[:]},nopython=True)")
        lines.append("def model_function(time, S, k):")
        lines.append("")
        lines.append("\tP = np.zeros(len(S))")
        lines.append("")
        lines.append("\t")
        lines.append("\t"+flow_profile_text)
        lines.append("")
        lines.append("")
        lines.append("\tidx = np.abs(F[0] - time).argmin()")
        lines.append("")
        lines.append("\tI = F[1:-1,idx]")
        lines.append("")
        lines.append('\tsigma_flow = F[-1,idx]')
        lines.append("")
        lines.append("")

        for m in model_text:
            lines.append(f"\t{m}")

        lines.append("\treturn P")
        lines.append("")
        lines.append("def wrapper_function(time, S, k):")
        lines.append("\treturn model_function(time, S, k)")
        lines.append("")

        lines.extend(self.write_variables_text())

        text = '\n'.join(lines)

        return text

    def write_model_matrix_text(self):
        '''

        Prototype for writing the model as an array.

        Parameters
        ----------

        Returns
        -------
        mat_text: str
            Rate equations in numpy matrix form.
        '''
        compounds = [*self.NetworkCompounds]
        reactions = [*self.NetworkReactions]

        species = self.species
        rate_consts = self.rate_constants
        inflows = self.inputs
        flow_ins = self.inflows
        flow_outs = self.outflows

        ratemat = [['0' for x in species] for x in species]

        for c in compounds:

            '''outgoing reactions'''
            for i in self.NetworkCompounds[c].In:

                reacs = self.NetworkReactions[i].Reactants
                ki = rate_consts[i]
                ind1 = compounds.index(c)
                if len(reacs) == 0:
                    token = f"(+{ki}*{inflows[c]})/{species[c]}"
                    ind2 = ind1

                if len(reacs) == 1:
                    n2 = reacs[0]
                    ind2 = compounds.index(n2)
                    token = '+' + ki

                if len(reacs) == 2:
                    n2 = reacs[0]
                    n3 = reacs[1]
                    ind2 = compounds.index(n3)
                    token = '+' + ki + '*' + species[n2]

                if len(reacs) == 3:
                    n2 = reacs[0]
                    n3 = reacs[1]
                    n4 = reacs[2]
                    ind2 = compounds.index(n3)
                    token = '+' + ki + '*' + species[n2] + '*' + species[n4]

                ratemat[ind1][ind2] += token

            for out in self.NetworkCompounds[c].Out:

                reacs = self.NetworkReactions[out].Reactants
                ind1 = compounds.index(c)
                ki = rate_consts[out]
                if len(reacs) == 1:
                    token = "-" + ki
                    ind2 = ind1

                if len(reacs) == 2:
                    z = reacs[:]
                    z.remove(c)
                    n2 = z[0]
                    ind2 = compounds.index(n2)
                    token = "-" + ki + '*' + species[c]


                if len(reacs) == 3:
                    z = reacs[:]
                    z.remove(c)
                    n2 = z[0]
                    n3 = z[1]
                    ind2 = compounds.index(n2)
                    token = "-" + ki + '*' + species[c] + '*' + species[n3]

                ratemat[ind1][ind2] += token

        mat_text = "["
        for r in ratemat:
            mat_text += "[" + ",".join(r) + "],\n"

        mat_text = mat_text.strip(",\n") + "]"

        return mat_text

def write_Jacobian_matrix_text(network):
    '''
    Prototype for writing the Jacobian matrix as text for the model.
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network to be written.

    Returns
    -------
    jac_text: str
        Jacobian matrix as text.
    '''
    compounds = [x for x in self.NetworkCompounds]
    reactions = [*self.NetworkReactions]

    species = self.species
    rate_consts = self.rate_constants
    inflows = self.inputs
    flow_ins = self.inflows
    flow_outs = self.outflows

    jac_mat = [['0' for x in species] for x in species]

    for c,comp1 in enumerate(compounds):
        for c2,comp2 in enumerate(compounds):
            element = ""
            for i in self.NetworkCompounds[comp1].In:
                if '_#0' in i:
                    pass
                elif comp2 in self.NetworkReactions[i].Reactants:
                    reacs = [species[x]
                                    for x in self.NetworkReactions[i].Reactants
                                        if x != comp2]
                    ki = f"+{rate_consts[i]}"
                    rctnt_elems = "*".join(reacs)
                    element += f"{ki}*{rctnt_elems}"
                else:
                    pass

            for o in self.NetworkCompounds[comp1].Out:
                if 'Sample' in o:
                    ki = f'-{flow_outs[o]}'
                    element += ki
                elif comp2 in self.NetworkReactions[o].Reactants:
                    reacs = [species[x]
                                    for x in self.NetworkReactions[o].Reactants
                                        if x != comp2]
                    ki = f"-{rate_consts[o]}"
                    rctnt_elems = "*".join(reacs)
                    element += f"{ki}*{rctnt_elems}"
                else:
                    pass

            jac_mat[c][c2] += element

    jac_text = ""
    for r in jac_mat:
        jac_text += "[" + ",".join(r) + "],\n"

    jac_text = jac_text.strip(",\n") + "]"

    return jac_text
