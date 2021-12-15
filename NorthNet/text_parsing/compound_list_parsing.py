from NorthNet import Classes
'''
For loading compound objects from a text file.
'''
def load_compounds_from_file(fname,
                            name_col = 0,
                            SMILES_col = 1,
                            delimiter = ','
                            ):
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
    with open(fname, "r") as file:
        for c,line in enumerate(file):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(delimiter)

                reagents[ ins[name_col] ] = Classes.Compound(ins[SMILES_col])

    return reagents
