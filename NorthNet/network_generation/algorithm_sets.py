from NorthNet import network_generation as n_gen
from NorthNet import Classes

def carbonyl_migration_isomers(network,
                               deprot_rule = Classes.Reaction_Template("deprotonation","[Ch:1][C:2]=[O:3].[O-:4] >> [C:1]=[C:2][O:3].[O-:4]",["[Ch:1][C:2]=[O:3","[O-:4]"],["[C:1]=[C:2][O:3]","[O-:4]"]),
                               prot_rule = Classes.Reaction_Template("protonation","[C:1]=[C:2][O:3].[O:4] >> [Ch:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"],["[Ch:1][C:2]=[O:3]","[O:4]"])):
    hydroxide = Classes.Compound("[OH-]")
    water = Classes.Compound("O")
    i = 0

    while i < 10:
        n_gen.extend_network_specific(network, [hydroxide], deprot_rule,[Classes.Substructure(x) for x in ["[C]=[C]([O])[O]","[C]=[C]=[C][O]","[C]=[C]=O","C=[C@@H]","C=[C@H]"]])

        n_gen.extend_network_specific(network, [water], prot_rule, [])

        n_gen.extend_network_specific(network, [hydroxide], deprot_rule,[Classes.Substructure(x) for x in ["[C]=[C]([O])[O]","[C]=[C]=[C][O]","[C]=[C]=O","C=[C@@H]","C=[C@H]"]])

        i += 1

def carbonyl_migration_isomers_stereo(network, deprot_rule = Classes.Reaction_Template("deprotonation","[Ch:1][C:2]=[O:3].[O-:4] >> [C:1]=[C:2][O:3].[O-:4]", ["[Ch:1][C:2]=[O:3]","[O-:4]"], ["[C:1]=[C:2][O:3]","[O-:4]"]),
                                        prot_rule_1 = Classes.Reaction_Template("protonation_@","[C:1]=[C:2][O:3].[O:4] >> [C@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@H:1][C:2]=[O:3]","[O:4]"]),
                                        prot_rule_2 = Classes.Reaction_Template("protonation_@@","[C:1]=[C:2][O:3].[O:4] >> [C@@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@@H:1][C:2]=[O:3]","[O:4]"])):
    hydroxide = Classes.Compound("[OH-]")
    water = Classes.Compound("O")

    i = 0
    while i < 20:
        n_gen.extend_network_specific(network, [hydroxide], deprot_rule,[])

        n_gen.extend_network_specific(network, [water], prot_rule_1, [])
        n_gen.extend_network_specific(network, [water], prot_rule_2, [])

        n_gen.extend_network_specific(network, [hydroxide], deprot_rule,[])

        i += 1

def carbonyl_migration_isomers_multiclass(network, deprot_rules = [Classes.Reaction_Template("deprotonation","[Ch:1][C:2]=[O:3].[O-:4] >> [C:1]=[C:2][O:3].[O-:4]", ["[Ch:1][C:2]=[O:3]","[O-:4]"], ["[C:1]=[C:2][O:3]","[O-:4]"])],
                                         prot_rules = [Classes.Reaction_Template("protonation_@","[C:1]=[C:2][O:3].[O:4] >> [C@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@H:1][C:2]=[O:3]","[O:4]"]),
                                         Classes.Reaction_Template("protonation_@@","[C:1]=[C:2][O:3].[O:4] >> [C@@H:1][C:2]=[O:3].[O:4]", ["[C:1]=[C:2][O:3]","[O:4]"], ["[C@@H:1][C:2]=[O:3]","[O:4]"])]):
    hydroxide = Classes.Compound("[OH-]")
    water = Classes.Compound("O")
    i = 0
    while i < 20:
        for d in deprot_rules:
            n_gen.extend_network_specific(network, [hydroxide], d, [])

        for p in prot_rules:
            n_gen.extend_network_specific(network, [water], p, [])

        i += 1
