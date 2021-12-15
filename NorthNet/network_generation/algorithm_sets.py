from NorthNet import Classes
from NorthNet import network_generation as n_gen

DEPROTONATION_RULE = "[Ch:1][C:2]=[O:3].[O-:4] >> [C:1]=[C:2][O:3].[O-:4]"
PROTONATION_RULE_0 ="[C:1]=[C:2][O:3].[O:4] >> [Ch:1][C:2]=[O:3].[O:4]"
PROTONATION_RULE_1 = "[C:1]=[C:2][O:3].[O:4] >> [C@H:1][C:2]=[O:3].[O:4]"
PROTONATION_RULE_2 = "[C:1]=[C:2][O:3].[O:4] >> [C@@H:1][C:2]=[O:3].[O:4]"

PROTONATION_REACTANT_SUBSTRUCTS = ["[C:1]=[C:2][O:3]","[O:4]"]

DEFAULT_DEPROTONATION_RULE = Classes.ReactionTemplate(
                                            "deprotonation",
                                            DEPROTONATION_RULE,
                                            ["[Ch:1][C:2]=[O:3","[O-:4]"],
                                            ["[C:1]=[C:2][O:3]","[O-:4]"]
                                            )
DEFAULT_PROTONATION_RULE = Classes.ReactionTemplate(
                                            "protonation",
                                            PROTONATION_RULE_0,
                                            PROTONATION_REACTANT_SUBSTRUCTS,
                                            ["[Ch:1][C:2]=[O:3]","[O:4]"]
                                            )

DEFAULT_PROTONATION_RULE_A = Classes.ReactionTemplate(
                                            "protonation_@",
                                            PROTONATION_RULE_1,
                                            ["[C:1]=[C:2][O:3]","[O:4]"],
                                            ["[C@H:1][C:2]=[O:3]","[O:4]"]
                                            )

DEFAULT_PROTONATION_RULE_B = Classes.ReactionTemplate(
                                            "protonation_@@",
                                            PROTONATION_RULE_2,
                                            ["[C:1]=[C:2][O:3]","[O:4]"],
                                            ["[C@@H:1][C:2]=[O:3]","[O:4]"]
                                            )

hydroxide = Classes.Compound("[OH-]")
water = Classes.Compound("O")

def carbonyl_migration_isomers(network,
                               deprot_rule = DEFAULT_DEPROTONATION_RULE,
                               prot_rule = DEFAULT_PROTONATION_RULE):

    i = 0
    while i < 10:
        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        n_gen.extend_network_specific(network, [water], prot_rule)

        n_gen.extend_network_specific(network, [hydroxide], deprot_rule)

        i += 1

def carbonyl_migration_isomers_stereo(
                                    network,
                                    deprot_rule = DEFAULT_DEPROTONATION_RULE,
                                    prot_rule_1 = DEFAULT_PROTONATION_RULE_A,
                                    prot_rule_2 = DEFAULT_PROTONATION_RULE_B):

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
                                                    DEFAULT_DEPROTONATION_RULE
                                                    ],

                                         prot_rules = [
                                                    DEFAULT_PROTONATION_RULE_A,
                                                    DEFAULT_PROTONATION_RULE_B
                                                    ]):
    i = 0
    while i < 20:
        for d_rule in deprot_rules:
            n_gen.extend_network_specific(network, [hydroxide], d_rule)

        for p_rule in prot_rules:
            n_gen.extend_network_specific(network, [water], p_rule)

        i += 1

