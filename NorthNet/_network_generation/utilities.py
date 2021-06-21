def cleanup(reactions):
    from rdkit import Chem
    
    '''Used to clear up the output of rdkit reaction function, parsing multiple
    reaction outcomes. Not perfect.'''

    reactions = [r.replace('[H+][O-]','O')    for r in reactions]
    reactions = [r.replace('[O-][H+]','O')    for r in reactions]
    #reactions = [r.replace('=[C@H]'   ,'=C')  for r in reactions]
    #reactions = [r.replace('=[C@@H]'   ,'=C') for r in reactions]
    reactions = [r.replace("/"   ,"") for r in reactions]
    reactions = [r.replace( r"\\"  ,"") for r in reactions]
    #reactions = [r.replace('@'   ,'') for r in reactions]

    reactions_out = []
    for c in range(0,len(reactions)):
        # Split the reaction up for inspection
        LHS = reactions[c].split(">>")[0].split(".")
        RHS = reactions[c].split(">>")[1].split(".")
        ins = ''
        for l in LHS:
            mol = Chem.MolFromSmiles(l)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        ins = ins.strip('.')
        ins += '>>'

        for r in RHS:
            mol = Chem.MolFromSmiles(r)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        reactions_out.append(ins.strip('.'))

    return reactions_out
