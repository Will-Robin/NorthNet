def name_abbreviation(name):
    '''
    Needs updating
    '''
    string_name = name[:]
    abbreviation = ''
    abbreviation += string_name[0]
    for x in string_name:
        if x == '_':
            abbreviation += string_name[string_name.index('_')+1]
            string_name =string_name[string_name.index('_')+1:]

    return abbreviations
