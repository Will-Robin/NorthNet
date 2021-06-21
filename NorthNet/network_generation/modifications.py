from NorthNet import Classes

def add_flow_terms(network, inputs):
    '''
    Add flow terms into network (inputs and outputs are
    empty species and )
    '''

    add_inputs = []
    for i in inputs:
        r_obj = Classes.ReactionInput("{}_#0>>{}".format(i,i))
        add_inputs.append(r_obj)

    network.add_inputs(add_inputs)

    add_outputs = []
    for c in network.NetworkCompounds:
        r_obj = Classes.ReactionOutput("{}>>Sample".format(c))
        add_outputs.append(r_obj)

    network.add_outputs(add_outputs)
