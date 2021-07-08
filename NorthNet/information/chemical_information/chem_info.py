'''
For importing stored chemical information on compounds.
'''

from pathlib import Path

path = Path(__file__)
script_dir = path.parent

# Load in file
info_container = []
with open(script_dir/'compound_properties.csv', 'r') as f:
    for c, line in enumerate(f):
        if c == 0:
            header = line.strip('\n').split(',')
        else:
            spl = line.strip('\n').split(',')
            info_container.append(spl)

# convert the lines into a dictionary
info_container = [list(i) for i in zip(*info_container)]
props_dict = {}
for n,i in zip(header, info_container):
    props_dict[n] = i

# create a dictionary assigning colours to compounds
colour_assignments = {k:v for k,v in
                              zip(props_dict['@ SMILES'], props_dict['colour'])}

for c,p in enumerate(props_dict['@@ SMILES']):
    if p in colour_assignments:
        pass
    else:
        colour_assignments[p] = props_dict['colour'][c]

# assign colours to other names
for c,p in enumerate(props_dict['Other_names']):
    for s in p.split(';'):
        if s == '':
            pass
        elif s in colour_assignments:
            pass
        else:
            colour_assignments[s] = props_dict['colour'][c]

# create a dict assgining names to smiles
canonical_SMILES = {k:v for k,v in
                       zip(props_dict['compound_name'], props_dict['@ SMILES'])}

for a,b in zip(props_dict['Other_names'], props_dict['@ SMILES']):
    for s in a.split(';'):
        if s != '':
            canonical_SMILES[s] = b

# for converting smiles to names
smiles_to_names = {}
for c,v in enumerate(props_dict['compound_name']):
    spl_name = v.split(' ')[0]
    smiles_to_names[props_dict['@ SMILES'][c]] = spl_name

# storing reaction class information
reaction_SMARTS = {}
reaction_class_colours = {}
reaction_class_short_names = {}
reaction_class_names = {}
with open(script_dir/'reaction_SMARTS_templates.tsv', 'r') as f:
    for c,line in enumerate(f):
        if c==0:
            pass
        else:
            ins = line.strip('\n').split('\t')
            reaction_SMARTS[ins[0]] = ins[3]
            reaction_class_colours[ins[0]] = ins[4]
            reaction_class_short_names[ins[0]] = ins[5]
            reaction_class_names[ins[0]] = ins[6]
