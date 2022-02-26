"""
Sketch of creating a text-based reaction score notation from Luca Cardelli (https://arxiv.org/abs/2005.08097)
"""

import math
from NorthNet import Classes

width = 10
horizontal_padding = 3
rxn_node_gap_pos = 1
up_ind = 0  # -1
down_ind = 2  # 1

spacer = " "
horizontal_line_token = "─"
reactant_down_token = "┯"
reactant_up_token = "┷"
product_pointer_up = "▲"
product_pointer_down = "▼"
vline_token = "│"


def substitute_string(string, new, position):
    """
    Replace the character at 'position' with 'new' in 'string'.

    Parameters
    ----------
    string: str

    new: str

    position: int

    Returns
    -------
    modified: str
    """

    modified = string[0:position] + new + string[position + 1 :]

    return modified


def create_reaction_block(network, reaction):

    compounds = [*network.NetworkCompounds]
    reactants = reaction.Reactants
    products = reaction.Products

    reaction_block = []
    for _ in compounds:
        reaction_block.append(
            [spacer * width, horizontal_line_token * width, spacer * width]
        )

    reactant_indices = []
    product_indices = []
    for reactant in reactants:
        index = compounds.index(reactant)
        reactant_indices.append(index)

    for product in products:
        index = compounds.index(product)
        product_indices.append(index)

    context_indices = reactant_indices + product_indices

    # find where to put the reaction node (a reference point)
    average_position = sum(context_indices) / len(context_indices)
    reaction_hpos = math.floor(width / 2) - 1
    rxn_vpos = math.floor(average_position)

    # draw lines for reactants
    for reac_count, r_i in enumerate(reactant_indices):

        vline_pos = reac_count * horizontal_padding + 1

        if r_i <= rxn_vpos:
            increment = 1
        else:
            increment = -1

        line_side = increment + 1
        for x in range(r_i + increment, rxn_vpos + increment, increment):
            for c, v in enumerate(reaction_block[x]):
                new_string = substitute_string(v, vline_token, vline_pos)
                reaction_block[x][c] = new_string

        if line_side == 0:
            line = reaction_block[x][0]
            newline = substitute_string(line, " ", vline_pos)
            reaction_block[x][0] = newline
            line = reaction_block[x][1]
            newline = substitute_string(line, "─", vline_pos)
            reaction_block[x][1] = newline
        if line_side == 2:
            line = reaction_block[x][-1]
            newline = substitute_string(line, " ", vline_pos)
            reaction_block[x][-1] = newline

        new_string = substitute_string(
            reaction_block[r_i][line_side], vline_token, vline_pos
        )
        reaction_block[r_i][line_side] = new_string

    # draw lines for products
    for prod_count, p_i in enumerate(product_indices):

        vline_pos = reaction_hpos + prod_count * horizontal_padding

        if p_i <= rxn_vpos:
            increment = 1
        else:
            increment = -1

        line_side = increment + 1
        for x in range(p_i + increment, rxn_vpos, increment):
            for c, v in enumerate(reaction_block[x]):
                new_string = substitute_string(v, vline_token, vline_pos)
                reaction_block[x][c] = new_string

        new_string = substitute_string(
            reaction_block[p_i][line_side], vline_token, vline_pos
        )
        reaction_block[p_i][line_side] = new_string

    # draw symbols for reactants
    for reac_count, r_i in enumerate(reactant_indices):
        vline_pos = reac_count * horizontal_padding + 1

        if r_i <= rxn_vpos:
            increment = 1
            r_token = reactant_down_token
            connect_token = "└"
            line_side = -1
        else:
            increment = -1
            r_token = reactant_up_token
            connect_token = "┌"
            line_side = -1

        new_line = substitute_string(
            reaction_block[rxn_vpos][line_side], connect_token, vline_pos
        )
        reaction_block[rxn_vpos][line_side] = new_line

        if vline_pos < reaction_hpos:
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos + 1,
            )
            reaction_block[rxn_vpos][line_side] = new_line
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos + 2,
            )
            reaction_block[rxn_vpos][line_side] = new_line
        else:
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos - 1,
            )
            reaction_block[rxn_vpos][line_side] = new_line
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos - 2,
            )
            reaction_block[rxn_vpos][line_side] = new_line

        compound_line = reaction_block[r_i][1]
        new_line = substitute_string(compound_line, r_token, vline_pos)
        reaction_block[r_i][1] = new_line

    # Draw symbol for products
    for prod_count, p_i in enumerate(product_indices):

        vline_pos = reaction_hpos + prod_count * horizontal_padding

        if p_i <= rxn_vpos:
            p_token = product_pointer_up
            p_token_pos = 2
            connect_token = "└"
        else:
            p_token = product_pointer_down
            p_token_pos = 0
            connect_token = "┐"

        line_side = -1
        new_line = substitute_string(
            reaction_block[rxn_vpos][line_side], connect_token, vline_pos
        )
        reaction_block[rxn_vpos][line_side] = new_line

        if vline_pos < reaction_hpos:
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos + 1,
            )
            reaction_block[rxn_vpos][line_side] = new_line
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos + 2,
            )
            reaction_block[rxn_vpos][line_side] = new_line
        else:
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos - 1,
            )
            reaction_block[rxn_vpos][line_side] = new_line
            new_line = substitute_string(
                reaction_block[rxn_vpos][line_side],
                horizontal_line_token,
                vline_pos - 2,
            )
            reaction_block[rxn_vpos][line_side] = new_line

        compound_line = reaction_block[p_i][p_token_pos]
        new_line = substitute_string(compound_line, p_token, vline_pos)
        reaction_block[p_i][p_token_pos] = new_line

    # draw reaction node
    insertion_line = reaction_block[rxn_vpos][2][:]
    new_line = substitute_string(insertion_line, "○", reaction_hpos)
    reaction_block[rxn_vpos][2] = new_line

    return reaction_block


def reaction_score_layout(network):
    """
    Create a 'reaction score' (https://arxiv.org/abs/2005.08097) from a
    network.

    Example:

    A + B → C

    B + C → A

    A + B → C + D

    A ──┯─────── ────────── ──┯───────
        │             ▲       │
        │             │       │
    B ──│──┯──── ──┯──│──── ──│──┯────
        └──○       └──○       └──○──┐
           ▼          │          ▼  │
    C ────────── ─────┷────  ───────│──
                                    │
                                    ▼
    C ────────── ──────────  ──────────

    Idea: each compound has a line (how to decide on ordering will become
    apparent. Each reaction has a horizontal block placed along the 'score'.
    Connections are made between reacting species as shown in the example.

    Let's try an ASCII diagram first.

    Parameters
    ----------

    Returns
    -------
    layout: ?

    """

    block_width = 10

    horizontal_padding = 2

    vertical_padding = 2

    node_column = 4
    node_row = 2
    reactant_column = 2
    product_column = 4

    compounds = [*network.NetworkCompounds]

    reactions = [*network.NetworkReactions]

    layout = []
    for c, reaction in enumerate(reactions):
        print(reaction)
        new_reaction_block = create_reaction_block(
            network, network.get_reaction(reaction)
        )
        layout.append(new_reaction_block)

    return layout


from NorthNet import text_parsing as file_import

reactions = """O=C(CO)CO.OC=C(O)CO>>O=C[C@@](O)(CO)C(O)(CO)CO
O=C(CO)CO.[OH-]>>OC=C(O)CO.[OH-]
O=C[C@H](O)CO.OC=C(O)CO>>O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO
O=C[C@H](O)CO.OC=C(O)CO>>O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO
O.OC=C(O)CO>>O.O=C[C@H](O)CO """

reaction_list = reactions.split("\n")

network = file_import.load_network_from_reaction_list(
    reaction_list, name="", description=""
)

layout = reaction_score_layout(network)
