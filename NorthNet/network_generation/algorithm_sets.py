from NorthNet import Classes
from NorthNet import network_generation as n_gen

default_deprotonation_rule = Classes.ReactionTemplate("deprotonation","[Ch:1][C:2]=[O:3].[O-:4] >> [C:1]=[C:2][O:3].[O-:4]",["[Ch:1][C:2]=[O:3","[O-:4]"],["[C:1]=[C:2][O:3]","[O-:4]"])
default_protonation_rule = Classes.ReactionTemplate("protonation","[C:1]=[C:2][O:3].[O:4] >> [Ch:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"],["[Ch:1][C:2]=[O:3]","[O:4]"])
default_protonation_rule_a = Classes.ReactionTemplate("protonation_@","[C:1]=[C:2][O:3].[O:4] >> [C@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@H:1][C:2]=[O:3]","[O:4]"])
default_protonation_rule_b = Classes.ReactionTemplate("protonation_@@","[C:1]=[C:2][O:3].[O:4] >> [C@@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@@H:1][C:2]=[O:3]","[O:4]"])

hydroxide = Classes.Compound("[OH-]")
water = Classes.Compound("O")

def carbonyl_migration_isomers(network,
                               deprot_rule = default_deprotonation_rule,
                               prot_rule = default_protonation_rule):

    i = 0

    while i < 10:
        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        n_gen.extend_network_specific(network, [water], prot_rule)

        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        i += 1

def carbonyl_migration_isomers_stereo(
                                    network,
                                    deprot_rule = default_deprotonation_rule,
                                    prot_rule_1 = default_protonation_rule_a,
                                    prot_rule_2 = default_protonation_rule_b):

    i = 0
    while i < 20:
        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        n_gen.extend_network_specific(network, [water], prot_rule_1)
        n_gen.extend_network_specific(network, [water], prot_rule_2)

        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        i += 1

def carbonyl_migration_isomers_multiclass(
                                        network,
                                        deprot_rules = [
                                                    default_deprotonation_rule
                                                    ],

                                         prot_rules = [
                                                    default_protonation_rule_a,
                                                    default_protonation_rule_b
                                                    ]):
    i = 0
    while i < 20:
        for d in deprot_rules:
            n_gen.extend_network_specific(network, [hydroxide], d)

        for p in prot_rules:
            n_gen.extend_network_specific(network, [water], p)

        i += 1
