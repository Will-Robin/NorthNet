"""
For loading reaction information fron text files.
"""
from NorthNet import Classes


def load_reaction_templates_from_file(fname, delimiter="\t"):
    """
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
    """

    lines = []
    with open(fname, "r") as file:
        for c, line in enumerate(file):
            lines = file.readlines()

    reaction_templates = {}
    for c, line in enumerate(lines):
        if c == 0:
            pass
        else:
            ins = line.strip("\n").split(delimiter)
            reaction_templates[ins[0]] = Classes.ReactionTemplate(
                ins[0], ins[3], ins[1].split("."), ins[2].split(".")
            )
    return reaction_templates
