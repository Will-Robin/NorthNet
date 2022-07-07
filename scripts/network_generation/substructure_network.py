"""
Script for establising Substructure Network functionality.
"""

from NorthNet import Classes
from NorthNet import network_visualisation as net_vis

reaction_file = "reactions_and_rules.tsv"

with open(reaction_file, "r") as file:
    text = file.read()

# parse text
lines = text.split("\n")

reaction_list = []
reaction_templates = []
seen_rules = []
for line in lines[1:]:  # skip header
    if line == "":
        continue
    entries = line.split("\t")

    reaction_smiles = entries[0]
    reaction_name = entries[1]
    reaction_smarts = entries[2]
    reactant_substructures = entries[3].split(".")
    product_substructures = entries[4].split(".")

    rule = Classes.ReactionTemplate(
        reaction_name, reaction_smarts, reactant_substructures, product_substructures
    )

    if rule.Name not in seen_rules:
        seen_rules.append(reaction_name)
        reaction_templates.append(rule)

    reaction_list.append(Classes.Reaction(reaction_smiles, reaction_template=rule))

# test adding reactions/creating substructure network from reactions
s_net = Classes.SubstructureNetwork(reaction_list, "", "")

# test removing compounds
compounds = [s_net.Compounds[c] for c in s_net.Compounds]
s_net.remove_compounds(compounds)

# test adding compounds
s_net.add_compounds(compounds)

# test creating substructure network from ReactionTemplates/ adding
# ReactionTemplates
s_net_2 = Classes.SubstructureNetwork([], "", "")
s_net_2.add_reaction_rules(reaction_templates)

# test removing reaction rules
s_net_2.remove_reaction_rules(reaction_templates)
