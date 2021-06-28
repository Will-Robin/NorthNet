def load_compounds_from_csv(fname, name_col = 0, SMILES_col = 1):
    '''
    Reads compounds from a .csv file.

    Parameters
    ----------
    fname: str or pathlib Path
        Path to the file containing nformation.
    name_col: int
        Column in the file which will give the keys for the output dict
    SMILES_col:
        Column containing the compound SMILES

    Returns
    -------
    reagents: dict
        Dictionary containing extracted compounds.
        {SMILES string: NorthNet Compound object}
    '''
    reagents = {}
    with open(fname, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(",")

                reagents[ ins[name_col] ] = Classes.Compound(ins[SMILES_col])

    return reagents

def load_reaction_templates_from_csv(fname, delimiter  = '\t'):
    '''
    Reads reaction templates from a .csv file.

    Assumed that the file is structure as follows:
    header\n
    name\treactant SMARTS\tproduct SMARTS\tReaction SMARTS
    etc.

    (ignores anything beyond 4th column)

    Parameters
    ----------
    fname: str
        file containing reaction templates and substructures.
    delimiter: str or pathlib Path
        Column delimiter for the file

    Returns
    -------
    reaction_templates: dict
        Dictionary of reaction templates.
        {reaction class name: NorthNet Reaction_Template}
    '''
    reaction_templates = {}
    with open(fname, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(delimiter)
                reaction_templates[ ins[0] ] = Classes.Reaction_Template(ins[0],
                                                            ins[3],
                                                            ins[1].split("."),
                                                            ins[2].split("."))

    return reaction_templates

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
