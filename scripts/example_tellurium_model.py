from NorthNet.Classes import Network
from NorthNet.Classes import ModelWriter
from NorthNet.text_parsing import load_network_from_reaction_list

reaction_file = "exampleReactionList.txt"

with open(reaction_file, "r") as file:
    text = file.read()

reaction_list = [x for x in text.split("\n") if x != ""]

network = load_network_from_reaction_list(reaction_list)

model = ModelWriter(network=network)

tellurium_model = model.to_tellurium_model()

print(tellurium_model)
