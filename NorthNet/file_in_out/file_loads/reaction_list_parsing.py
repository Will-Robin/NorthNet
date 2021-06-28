def load_compounds_from_csv(fname):
    '''
    Reads compounds from a .csv file.

    Parameters
    ----------
    fname: str or pathlib Path
        Path to the file containing nformation.
    Returns
    -------
    reagents: dict
        Dictionary containing extracted compounds. SMILES: NorthNet Compound object
    '''
    reagents = {}
    with open(fname, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(",")

                reagents[ ins[0] ] = Classes.Compound(ins[1])

    return reagents

def get_reaction_templates(fname, delimiter  = '\t'):
    '''
    Reads reaction templates from a .csv file.

    Parameters
    ----------
    fname: str
        file containing reaction templates and substructures.
    Returns
    -------
    reaction_templates: dict
        Dictionary of reaction templates. reaction class: NorthNet Reaction_Template
    '''
    reaction_templates = {}
    with open(fname, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(delimiter)
                reaction_templates[ ins[0] ] = Classes.Reaction_Template(ins[0], ins[3], ins[1].split("."), ins[2].split("."))

    return reaction_templates

def get_reaction_list_from_database(file, limit = 2e200):
    '''
    Creates reaction_list of Reaction objects from a data file.
    The limit is in place so that the whole file does not have to be imported
    Parameters
    ----------
    file: str
        Path to the file with reaction information.
    limit: int
        Maximum number of lines to be read.

    Returns
    -------
    reaction_list: list
        A list of NorthNet_functions Reaction objects derived from file.
    '''
    arr = []
    n=0
    with open(file,'r') as f:
        for line in f:
            z = line.strip('\n')
            arr.append(z.split('\t'))
            n+=1
            if n == limit:
                break

    header = arr[0]
    data = arr[1:]

    database_entries = []
    for d in data:
        x = Classes.Reaction_Database_Entry(header,d)
        database_entries.append(x)

    return header, database_entries

def read_mapped_reactions_file(fname):
    mapped_reactions = {}

    with open(fname, 'r') as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip('\n').split(',')
                mapped_reactions[ins[0]] = ins[1]

    return mapped_reactions

def load_schemes(scheme_dir):
    home = os.getcwd()
    os.chdir(scheme_dir)
    # get reaction lists
    schemes = {}
    for file in os.listdir():
        if file.endswith(".txt"):
            reactions = []
            with open(file, "r") as f:
                for line in f:
                    reactions.append(line.strip("\n"))
            schemes[file.split("_")[0]] = reactions
    os.chdir(home)
    return schemes

def get_reaction_clusters(fname):
    rxn_clusters = {}
    process_line = lambda x: x.strip("\n").split(",")
    with open(fname, "r") as f:
        for line in f:
            ins = process_line(line)
            rxn_clusters[ins[0]] = [x for x in ins[1:] if x != ""]
    return rxn_clusters

def reactions_from_list(reaction_list):

    header = ["Reaction", "Description", "Reaction ID","References"]
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    ID_dict = {S:str(c) for c,S in enumerate(alphabet)}

    add_to_net = []
    for r in reaction_list:

        reactants = r.split(">>")[0].split(".")
        products = r.split(">>")[1].split(".")

        r_mols = [Chem.MolFromSmiles(x) for x in reactants]
        p_mols = [Chem.MolFromSmiles(x) for x in products]

        r_canon = [Chem.MolToSmiles(x, isomericSmiles = True) for x in r_mols]
        p_canon = [Chem.MolToSmiles(x, isomericSmiles = True) for x in p_mols]

        r_string = r_ops.reconstitute_reaction(r_canon,p_canon)

        rID = "".join([ID_dict[s] for s in r_string if s in ID_dict])

        line = [r_string, "User input reaction.", rID, "None"]

        add_to_net.append(Classes.Reaction( Classes.Reaction_Database_Entry(header, line) ))

    return add_to_net
