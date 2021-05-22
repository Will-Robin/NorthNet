from pathlib import Path

def readTwoColInfo(file, col_idx_a, col_idx_b):
    assignments = {}
    with open(file, 'r') as f:
        for c, line in enumerate(f):
            if c > 0:
                spl = line.strip('\n').split(',')
                assignments[spl[col_idx_a]] = spl[col_idx_b]
            else:
                pass
    return assignments

path = Path(__file__)
script_dir = path.parent

info_container = []
with open(script_dir/'compound_properties.csv', 'r') as f:
    for c, line in enumerate(f):
        if c == 0:
            header = line.strip('\n').split(',')
        else:
            spl = line.strip('\n').split(',')
            info_container.append(spl)

info_container = [list(i) for i in zip(*info_container)]

props_dict = {}
for n,i in zip(header, info_container):
    props_dict[n] = i

colour_assignments = {k:v for k,v in
                              zip(props_dict['@ SMILES'], props_dict['colour'])}

for c,p in enumerate(props_dict['@@ SMILES']):
    if p in colour_assignments:
        pass
    else:
        colour_assignments[p] = props_dict['colour'][c]

molecular_masses = {k:float(v) for k,v in
                           zip(props_dict['@ SMILES'], props_dict['Mr_gmol-1'])}
canonical_SMILES = {k:v for k,v in
                       zip(props_dict['compound_name'], props_dict['@ SMILES'])}

for a,b in zip(props_dict['Other_names'], props_dict['@ SMILES']):
    for s in a.split(';'):
        if s != '':
            canonical_SMILES[s] = b

smiles_to_names = {}
init_can_SMILES = [c for c in canonical_SMILES]
for c in init_can_SMILES:
    spl_name = c.split(' ')
    canonical_SMILES[spl_name[0]] = canonical_SMILES[c]
    smiles_to_names[canonical_SMILES[c]] = spl_name[0]

class_assignments =  {k:v for k,v in
                               zip(props_dict['@ SMILES'], props_dict['Class'])}
for sm,cls in zip(props_dict['@@ SMILES'], props_dict['Class']):
    class_assignments[sm] = cls

reaction_colours = readTwoColInfo(script_dir/'reaction_colour_assignments.csv', 0, 1)
fr_assign = readTwoColInfo(script_dir/'fragment_assignments.csv', 0, 1)
frag_assignments = {float(f):fr_assign[f] for f in fr_assign}
frag_colours = readTwoColInfo(script_dir/'fragment_assignments.csv', 0, 2)
