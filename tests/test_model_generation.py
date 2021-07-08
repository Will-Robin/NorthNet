import os
import numpy as np
from NorthNet import Classes
import pickle

# load the test data report
data_report = Classes.DataReport(file = 'tests/Test_data_report.csv')
# load test reactions
with open('tests/Test_reaction_list.txt', 'r') as f:
    lines = f.readlines()
reactions = [l.strip('\n') for l in lines]

with open("tests/FormoseReactionNetwork.pickle", 'rb') as f:
    formose = pickle.load(f)

model_reactions = [formose.NetworkReactions[r] for r in reactions]

model_network = Classes.Network(model_reactions,
                            'test_model',
                            'for testing model writing'
)

model = Classes.ModelWriter(
                        network = model_network,
                        experiment = data_report,
                        input_token = '_#0',
                        output_token = 'Sample',
                        flowrate_time_conversion = 3600,
                        time_limit = False,
                        lead_time = data_report.series_values[0]-1000)

model_text = model.write_to_module_text(numba_decoration = True)

output_file = 'tests/test_model.py'
with open(output_file,'w') as f: f.write(model_text)
