from rdkit import Chem
from rdkit.Chem import AllChem
from NorthNet import Classes
import NorthNet.molecule_operations as mol_ops


def run_rdkit_reaction(reactant_compounds, reaction_template):
    """
    Performs a chemical reaction.

    Mapping can be added:
        Set atom map numbers of reactants before performing the reaction.
        Atom mappings can be removed:
        rdkit.Chem.rdChemReactions.RemoveMappingNumbersFromReactions(
                                                 (ChemicalReaction)reaction
                                                 ) → None

    Parameters
    ----------
    reactant_compounds: tuple
        tuple of NorthNet Compound objects which take part in the reaction.
    reaction_template: list
        NorthNet Reaction_Template object.

    Returns
    -------
    reactions: list
        A list of NorthNet Reaction objects.
    """

    reactions = []
    reactants = tuple(r.Mol for r in reactant_compounds)

    # Run the reaction to give a list of products sets
    product_sets = reaction_template.Reaction.RunReactants(reactants)

    for product_set in product_sets:
        for product in product_set:
            # The products are checked for chiral information which has been
            # transferred to achiral carbons, which is removed.
            product = mol_ops.incorrect_chiral_H_solve(product)
            Chem.SanitizeMol(product)

        rxn = AllChem.ChemicalReaction()  # Create an empty chemical reaction
        [rxn.AddReactantTemplate(r) for r in reactants]
        [rxn.AddProductTemplate(p) for p in product_set]

        reaction_SMILES = AllChem.ReactionToSmiles(rxn)

        reactions.append(
            Classes.Reaction(reaction_SMILES, reaction_template=reaction_template)
        )

    return reactions
