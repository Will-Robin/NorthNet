def convert_dataset_to_vector(dataset, header):

    matrix = np.zeros((len(dataset.time),len(header)))
    for d in dataset.dependents:
        ind = header.index(d)
        for x in range(0,len(dataset.time)):
            matrix[x,ind] = dataset.dependents[d][x]

    vectors = {}
    for x in range(0,len(matrix)):
        vectors[x] = matrix[x]

    return vectors

def convert_data_to_vector(vector_header, data):
    vec = np.zeros(len(vector_header))
    for d in data:
        if d == "C=O M" or d == "O=C(CO)CO M":
            pass
        else:
            vec[vector_header.index(d)] = data[d]
    return vec

def convert_dataset_to_vector(dataset, header):

    matrix = np.zeros((len(dataset.time),len(header)))
    for d in dataset.dependents:
        ind = header.index(d)
        for x in range(0,len(dataset.time)):
            matrix[x,ind] = dataset.dependents[d][x]

    vectors = {}
    for x in range(0,len(matrix)):
        vectors[x] = matrix[x]

    return vectors

def matrix_from_data_set_dict(data_sets,header):
    from NetFit import fitting_functions as ff

    pos = 0
    all_points = 0
    for d in data_sets:
        all_points += len(data_sets[d].time)

    data_mat = np.zeros((all_points,len(header)))
    labels = []
    for c,d in enumerate(data_sets):
        vec = ff.convert_dataset_to_vector(data_sets[d],header)
        for c2,v in enumerate(vec):
            data_mat[pos+c2] = vec[v]
            labels.append("{};{}".format(d,data_sets[d].time[c2]))
        pos += len(data_sets[d].time)
    return data_mat,labels
