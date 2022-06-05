from NorthNet.Classes import Network
from NorthNet.Classes import ModelWriter
from NorthNet.Classes import DataReport
from NorthNet.Loading import load_network_from_reaction_list

reaction_file = "exampleReactionList.txt"
data_file = "exampleData.csv"

with open(reaction_file, "r") as file:
    text = file.read()

reaction_list = [x for x in text.split("\n") if x != ""]

network = load_network_from_reaction_list(reaction_list)
experiment = DataReport(file=data_file)

model = ModelWriter(network=network, experiment=experiment)

tellurium_model = model.to_tellurium_model()

print(tellurium_model)
