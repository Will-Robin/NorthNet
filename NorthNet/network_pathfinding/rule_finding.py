def create_rule_subunit(rules):
    subunit = nx.DiGraph()
    subunit.add_node('input', name = 'input')
    subunit.add_node('output', name = 'output')
    for d in rules:
        units = d.split('-')
        at_counter = 1
        for c,x in enumerate(units):

            if '@' in x:
                node = '{}{}{}'.format(x,d,at_counter)
                at_counter += 1
            else:
                node = '{}'.format(x)

            subunit.add_node(node, name = x)

            if c == 0:
                subunit.add_edge('input', node)

            if c == len(units)-1:
                subunit.add_edge(node,'output')
            else:
                if '@' in units[c+1]:
                    dest_node = '{}{}{}'.format(units[c+1],d,at_counter)
                else:
                    dest_node = '{}'.format(units[c+1])

                subunit.add_edge(node,dest_node)

    return subunit

def create_subunit_path(rule):
    path = nx.DiGraph()
    path.add_node('input', name = 'input')
    path.add_node('output', name = 'output')

    units = rule.split('-')
    at_counter = 1
    for c,x in enumerate(units):

        if '@' in x:
            node = '{}{}{}'.format(x,rule,at_counter)
            at_counter += 1
        else:
            node = '{}'.format(x)

        path.add_node(node, name = x)

        if c == 0:
            path.add_edge('input', node)

        if c == len(units)-1:
            path.add_edge(node,'output')
        else:
            if '@' in units[c+1]:
                dest_node = '{}{}{}'.format(units[c+1],rule,at_counter)
            else:
                dest_node = '{}'.format(units[c+1])

            path.add_edge(node,dest_node)

    return path
