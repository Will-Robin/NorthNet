import yaml
from rdkit import Chem

from NorthNet import Classes
from NorthNet import text_parsing as text_p
from NorthNet import network_generation as n_gen


def generate_epimers(network, deprotonation_rules=[], protonation_rules=[]):
    """
    A function to generate all of the epimers of a sugar: over-iterates
    multiple times through a series of protonation/deprotonation reactions.

    Parameters
    ----------
    network: NorthNet.Network

    deprotonation_rules: list[NorthNet.ReactionTemplate]

    protonation_rules: list[NorthNet.ReactionTemplate]

    Returns
    -------
    None
    """

    hydroxide = Classes.Compound("[OH-]")
    water = Classes.Compound("O")

    i = 0
    reaction_number = len(network.NetworkReactions)
    while i < 0:
        for d_rule in deprotonation_rules:
            n_gen.extend_network_specific(network, [hydroxide], d_rule)

        for p_rule in protonation_rules:
            n_gen.extend_network_specific(network, [water], p_rule)

        new_reaction_number = len(network.NetworkReactions)

        i = new_reaction_number - reaction_number

        reaction_number = new_reaction_number

# Load in some variables
with open('params.yaml', 'r') as file:
    text = file.read()

info = yaml.load(text, Loader = yaml.FullLoader)

"""Get reaction components"""
# there may be an error such as:
# 'mapped atoms in the reactants were not mapped in the products.
# unmapped numbers are: 5'
# the error is attributable to the Cannizzaro reaction SMARTS, in which
# O:5 is not mapped to the products

reaction_SMARTS_file = info["reaction-smarts-file"]

reactions = text_p.load_reaction_templates_from_file(reaction_SMARTS_file)

C_patt = Chem.MolFromSmarts("[C]")
count_carbons = lambda x: x.GetSubstructMatches(C_patt)

"""Name"""
network_name = info["network-name"]
description = info["network-description"]

"""Boundary conditions"""
# iterations  overshoot for C6, but do so to get
# all reaction paths and compounds possible up
# to C6 compounds
iterations = info["iterations"]
start_smiles = info["initiator-smiles"]

initiator_species = [Classes.Compound(x) for x in start_smiles]
reaction_network = Classes.Network([], network_name, description)

reaction_network.add_compounds(initiator_species)

"""Reactivity constructor"""
reaction_pattern = info["reaction-rules"]
deprotonation_rules = [x for x in reaction_pattern if "deprotonation" in x]
protonation_rules = [x for x in reaction_pattern if "protonation" in x]

"""Expansion operation"""
x = 0
while x < iterations:
    for task in reaction_pattern:
        n_gen.extend_network_task(reaction_network, reactions[task])

    generate_epimers(
        reaction_network,
        deprotonation_rules=[reactions[d] for d in deprotonation_rules],
        protonation_rules=[reactions[p] for p in protonation_rules],
    )

    """Removing compounds > C6"""
    # In effect, this is equivalent to setting all chain-growing reaction
    # rules to not occur for C6 compounds
    # i.e. [$(C(O)=CO)!$(C(O)=C(O)C(O)C(O)C(O)CO)], etc.
    remove_compounds = [
        reaction_network.NetworkCompounds[c]
        for c in reaction_network.NetworkCompounds
        if len(count_carbons(reaction_network.NetworkCompounds[c].Mol)) > 6
    ]
    reaction_network.remove_compounds(remove_compounds)

    x += 1

compound_number = len(reaction_network.NetworkCompounds)
reaction_number = len(reaction_network.NetworkReactions)

print(f"Generated {compound_number} compounds and {reaction_number} reactions.")
