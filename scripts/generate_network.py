import sys
from rdkit import Chem
from NorthNet import Classes
from NorthNet import text_parsing as text_p
from NorthNet import network_generation as n_gen

def generate_epimers(
                    network,
                    deprotonation_rules = [],
                    protonation_rules = []
                    ):
    i = 0
    reaction_number = len(network.NetworkReactions)
    while i < 0:
        for d_rule in deprot_rules:
            n_gen.extend_network_specific(network, [hydroxide], d_rule)

        for p_rule in prot_rules:
            n_gen.extend_network_specific(network, [water], p_rule)

        new_reaction_number = len(network.NetworkReactions)

        i = new_reaction_number - reaction_number

        reaction_number = new_reaction_number


'''Get reaction components'''
# there may be an error such as:
# 'mapped atoms in the reactants were not mapped in the products.
# unmapped numbers are: 5'
# the error is attributable to the Cannizzaro reaction SMARTS, in which
# O:5 is not mapped to the products

reaction_SMARTS_file = 'reaction_SMARTS_templates.tsv'

reactions = text_p.load_reaction_templates_from_file(reaction_SMARTS_file)

C_patt = Chem.MolFromSmarts("[C]")
count_carbons = lambda x: x.GetSubstructMatches(C_patt)

'''Name'''
network_name = "formose-reaction"
description = "A formose reaction generated up to C6 compounds."

'''Boundary conditions'''
# iterations  overshoot for C6, but do so to get
# all reaction paths and compounds possible up
# to C6 compounds
iterations = 6
start_smiles = ['O','[OH-]','O=CCO', 'C=O', 'OC=CO','OC1COC(CO1)O']

initiator_species = [Classes.Compound(x) for x in start_smiles]
reaction_network = Classes.Network([], network_name, description)

reaction_network.add_compounds(initiator_species)

'''Reactivity constructor'''
# Choose reaction rules to include in the network generation algorithm.
deprotonation_rules = ['deprotonation-t1','deprotonation-t2',
                       'deprotonation-t3','deprotonation-t4',
                       'deprotonation-t5','deprotonation-t6',
                       'deprotonation-t7','deprotonation-t8']

protonation_rules = ['protonation-t1-@','protonation-t1-@@',
                     'protonation-t2-@','protonation-t2-@@',
                     'protonation-t3-@','protonation-t3-@@',
                     'protonation-t4-@','protonation-t4-@@',
                     'protonation-t5-@','protonation-t5-@@',
                     'protonation-t6-@','protonation-t6-@@',
                     'protonation-t7-@','protonation-t7-@@',
                     'protonation-t8-@','protonation-t8-@@']

aldol_addition_rules = ['aldol_addition_fast_pt_C=O-t1-@',
                        'aldol_addition_fast_pt_C=O-t1-@@',
                        'aldol_addition_fast_pt_C=O-t2-@',
                        'aldol_addition_fast_pt_C=O-t2-@@',
                        'aldol_addition_fast_pt_C=O-t3-@',
                        'aldol_addition_fast_pt_C=O-t3-@@',
                        'aldol_addition_fast_pt_C=O-t4-@',
                        'aldol_addition_fast_pt_C=O-t4-@@',
                        'aldol_addition_fast_pt_C=O-t5-@',
                        'aldol_addition_fast_pt_C=O-t5-@@',
                        'aldol_addition_fast_pt_C=O-t6-@',
                        'aldol_addition_fast_pt_C=O-t6-@@',
                        'aldol_addition_fast_pt_C=O-t7-@',
                        'aldol_addition_fast_pt_C=O-t7-@@',
                        'aldol_addition_fast_pt_C=O-t8-@',
                        'aldol_addition_fast_pt_C=O-t8-@@']

sugar_aldol_rules = ['aldol_sugar_addition_fast_pt-t1-@@-@',
                     'aldol_sugar_addition_fast_pt-t2-@@-@',
                     'aldol_sugar_addition_fast_pt-t1-@@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@-@',
                     'aldol_sugar_addition_fast_pt-t1-@-@',
                     'aldol_sugar_addition_fast_pt-t1-@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@-@@']

others = ['retroaldol_to_enol-t1','retroaldol_to_enol-t2',
            'Cannizzaro', 'dimer_dissociation']

reaction_pattern = deprotonation_rules + protonation_rules + \
                    aldol_addition_rules + sugar_aldol_rules + \
                    others

'''Expansion operation'''
x = 0
while x < iterations:
    for task in reaction_pattern:
        n_gen.extend_network_task(reaction_network, reactions[task])

    generate_epimers(
                    reaction_network,
                    deprotonation_rules = [reactions[d] for d in deprotonation_rules],
                    protonation_rules = [reactions[p] for p in protonation_rules]
                    )

    '''Removing compounds > C6'''
    # In effect, this is equivalent to setting all chain-growing reaction
    # rules to not occur for C6 compounds
    # i.e. [$(C(O)=CO)!$(C(O)=C(O)C(O)C(O)C(O)CO)], etc.
    remove_compounds = [reaction_network.NetworkCompounds[c]
                                    for c in reaction_network.NetworkCompounds
        if len(count_carbons(reaction_network.NetworkCompounds[c].Mol)) > 6]
    reaction_network.remove_compounds(remove_compounds)

    x+=1

compound_number = len(reaction_network.NetworkCompounds)
reaction_number = len(reaction_network.NetworkReactions)

print(f'Generated {compound_number} compounds and {reaction_number} reactions.')


