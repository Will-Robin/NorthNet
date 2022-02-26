import json


def dictionary_from_layout(json_string):
    """

    Parameters
    ----------
    json_string: str
        JSON string from dot layout.
    Returns
    -------
    pos: dict
        Dictionary of x,y coordinates.
    """
    layout_in_string = json.loads(json_string)

    pos = {}
    for l_obj in layout_in_string["objects"]:
        pos[l_obj["name"]] = [float(x) for x in l_obj["pos"].split(",")]

    return pos
