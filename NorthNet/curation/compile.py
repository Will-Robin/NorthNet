def compile_database(in_directory = r"C:\Users\willi\Documents\PROJECTS\NPC_Project\Databases\Reaxys_Reaction_Database\Individual_files", out_directory = r'C:\Users\willi\Documents\PROJECTS\NPC_Project\Databases\Reaxys_Reaction_Database\Compiled_Database'):
    '''
    Parameters
    ----------
    in_directory: str
        path to directory holding individual data files.
    out_directory: str
        path to directory in which to place complie database
    Returns
    -------
    None
    '''
    timestr = time.strftime("%Y%m%d-%H%M%S")

    os.chdir(in_directory)
    list = os.listdir()

    ALL = []
    with open(list[10],'r') as f:
        ALL.append(f.readline())
        pass
    print(ALL)
    for file in list:
        if file.endswith('.xls'):
            print(file)
            with open(file,'r') as f:
                next(f)
                for line in f:
                    ALL.append(line)

    os.chdir(out_directory)

    with open('{}_All_reactions_reaxys.txt'.format(timestr), 'w') as out:
        for a in ALL:
            out.write(a)
