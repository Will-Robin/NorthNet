from NorthNet.Classes import DataReport
from NorthNet.Classes import ReactionInput
from NorthNet.Classes import ReactionOutput
from NorthNet.Classes import Network
from NorthNet.Classes import ModelWriter
from NorthNet.text_parsing import load_network_from_reaction_list

data = DataReport(file = "exampleData.csv")

reaction_file = "exampleReactionList.txt"

with open(reaction_file, "r") as file:
    text = file.read()

reaction_list = [x for x in text.split("\n") if x != ""]

network = load_network_from_reaction_list(reaction_list)

for comp in network.NetworkCompounds:
    output = ReactionOutput(f"{comp}>>#0")
    network.add_output_process(output)

input_compounds = ["O=C(CO)CO", "[OH-]", "C=O"]

for smiles in input_compounds:
    input = ReactionInput(f"{smiles}_#0>>{smiles}")
    network.add_input_process(input)

model = ModelWriter(network=network, experiment = data, lead_time = data.series_values[0])

model_text = model.write_to_module_text(numba_decoration="compile")

with open("exampleModel.py", "w") as file:
    file.write(model_text)

