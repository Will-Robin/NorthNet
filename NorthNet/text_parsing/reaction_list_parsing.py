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
    from NorthNet import Classes

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

    from NorthNet import Classes

    with open(fname, "r") as f:
        for c,line in enumerate(f):
            lines = f.readlines()
            
    reaction_templates = {}
    for c,line in enumerate(lines):
        if c == 0:
            pass
        else:
            ins = line.strip("\n").split(delimiter)
            reaction_templates[ ins[0] ] = Classes.Reaction_Template(ins[0],
                                                        ins[3],
                                                        ins[1].split("."),
                                                        ins[2].split("."))
    return reaction_templates
