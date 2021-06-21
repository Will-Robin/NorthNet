def import_file_section(file, start_token, end_token):
    '''
    Parameters
    ----------
    file: path to file
    start_token: str
        String in line to start reading file from.
    end_token:
        String in line to end reading file from.
    '''

    spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
    readstate = False
    c_set = []
    with open(file, 'r', encoding = 'latin1') as f:
        for c,line in enumerate(f):
            if start_token in line:
                readstate = True
                line = next(f)
            if end_token in line:
                readstate = False
            if readstate:
                newline = spl_lin(line)
                c_set.append(newline)

    return c_set
