def create_image(smiles, name, fig_size = (10,10)):
    '''
    Creates an image from a SMILES structure and saves it to the directory.
    Parameters
    ----------
    smiles: str
        SMILES structure of the molecule to be plotted.
    name: str
        name for the output file.
    Returns
    -------
    fname: str
        name of the written file.
    '''
    fname = '{}.png'.format(name).replace("@","").replace(":","").replace("[","").replace("]","")

    m = Chem.MolFromSmiles(smiles,sanitize = True)

    if m == None:
        print("Could not parse as SMILES, attempting to parse as SMARTS.")
        m = Chem.MolFromSmarts(smiles)
    else:
        pass

    for atom in m.GetAtoms():
        atom.SetAtomMapNum(0)

    AllChem.Compute2DCoords(m)

    drawer = rdMolDraw2D.MolDraw2DSVG(fig_size[0],fig_size[1])

    opts = drawer.drawOptions()

    DrawingOptions.atomLabelFontSize = 55
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 3.0
    DrawingOptions.atomLabelFontFace = 'Arial'
    opts.clearBackground=False

    Draw.MolToFile( m, fname )
    Draw.MolToFile( m, "temp.svg" )
    cairosvg.svg2png( url='./temp.svg', write_to= fname )

    return "{}".format(fname)

def draw_image_2(mol, substruct, fname = 'test.png'):
    from rdkit.Chem import rdDepictor
    print('mol',Chem.MolToSmiles(mol), 'substr',Chem.MolToSmarts(substruct))
    Chem.SanitizeMol(mol)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
    drawer.DrawMolecule(mol,highlightAtoms=mol.GetSubstructMatch(substruct))
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    with open('temp.svg', 'w') as f:
        f.write(svg)
    cairosvg.svg2png( url='temp.svg', write_to= fname )
