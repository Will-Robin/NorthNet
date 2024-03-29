from NorthNet import Classes
from NorthNet.reaction_operations import conversions


def load_network_from_reaction_list(reaction_list, name="", description=""):
    """
    Create a NorthNet Network from a list of reactions

    Parameters
    ----------
    reaction_list: list[str]
        list of reaction SMILES strings
        Format: e.g. C=O.OC=C(O)CO>>O=C([C@@H](CO)O)CO

    Returns
    -------
    network: NorthNet.Classes.Network
    """

    rxns = []
    for reaction in reaction_list:
        rxns.append(Classes.Reaction(reaction))

    network = Classes.Network(rxns, name, description)

    return network


def node_edge_list_from_gdf(filename):
    """
    Get network edges and coordinates from .gdf file

    Parameters
    ----------
    file: str or pathlib Path
        Path to file

    Returns
    -------
    node_list, edge_list: list[dict]

    """
    with open(filename, "r") as file:
        coordinates = []
        node_list_header = []

        for c, line in enumerate(file):
            if c == 0:
                spl = line.strip("\n").split(",")
                node_list_header = [x.split(" ")[-2] for x in spl]
                continue

            if "edgedef" in line:
                break

            spl = line.strip("\n").split(",")
            coordinates.append({k: v for k, v in zip(node_list_header, spl)})

    node_list = {
        k["name"].strip(" "): [float(k["x"]), float(k["y"])] for k in coordinates
    }

    with open(filename, "r") as file:
        edge_list = []
        edge_list_header = []

        readstate = False
        for c, line in enumerate(file):
            if "edgedef" in line:
                spl = line.strip("\n").strip(">edgedef").split(",")
                edge_list_header = spl

            if readstate:
                spl = line.strip("\n").split(",")
                edge_list.append({k: v for k, v in zip(edge_list_header, spl)})

    return node_list, edge_list
