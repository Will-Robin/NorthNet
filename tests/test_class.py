from NorthNet import Classes
from rdkit.Chem import AllChem

rxn = AllChem.ChemicalReaction()
rt = Classes.ReactionTemplate('C=O', 'C=O>>CO', 'C=O', 'CO')

Classes.Reaction(rxn, rt, info = {})

print('''class ReactionInput:
        reaction_input_string arg must be str like SMILES_#0.''')
