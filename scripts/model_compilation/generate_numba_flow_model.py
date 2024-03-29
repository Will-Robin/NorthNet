import numpy as np
from NorthNet.Classes import DataReport
from NorthNet.Classes import InputProcess
from NorthNet.Classes import OutputProcess
from NorthNet.Classes import Network
from NorthNet.Classes import ModelWriter
from NorthNet.Loading import load_network_from_reaction_list

data = DataReport(file="exampleData.csv")
# Convert units into SI base units, given knowledge of the units in the file.
data.conditions["reactor_volume/ uL"] /= 1e6
for cond in data.conditions:
    if "flow" in cond and not "time" in cond:
        data.conditions[cond] *= 1e-6/3600

reaction_file = "exampleReactionList.txt"

with open(reaction_file, "r") as file:
    text = file.read()

reaction_list = [x for x in text.split("\n") if x != ""]

network = load_network_from_reaction_list(reaction_list)

model = ModelWriter(network=network, experiment=data)

# Trim the flow profiles down to limit the amount of compiling needed
min_time = data.series_values[0] - 1000
max_time = data.series_values[-1]
inds = np.where(
    (model.flow_profile_time > min_time) & (model.flow_profile_time < max_time)
)[0]

model.flow_profile_time = model.flow_profile_time[inds]
model.flow_profile_time -= model.flow_profile_time[0]
model.sigma_flow = model.sigma_flow[inds]
for fl in model.flow_profiles:
    model.flow_profiles[fl] = model.flow_profiles[fl][inds]

model_text = model.write_to_module_text(numba_decoration="compile")

print(model_text)
