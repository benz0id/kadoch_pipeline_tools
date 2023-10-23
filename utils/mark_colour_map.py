"""
cols = [
    "orange",
    "limegreen",
    "mediumblue",
    "fuchsia",
    "chocolate",
    "darkgreen",
    "slateblue",
    "mediumpurple",
    "grey",
    "grey",
    "deepskyblue",
    "darkmagenta",
    "teal"
]
"""
from typing import Dict

DEFAULT_COLOUR_MAP = {'9ala': 'green',
                      'arid1a': 'orange',
                      'arid1b': '"#FF8C00"',
                      'arid2': 'gold',
                      'atac': 'red',
                      'bicra': '"#00BF86"',
                      'bicral': '"#00BF86"',
                      'brd4': 'green',
                      'brd9': 'green',
                      'brg1': 'blue',
                      'c1': '"#0091FF"',
                      'ctcf': '"#BF49A6"',
                      'empty': 'grey',
                      'gltscr': 'green',
                      'gltscr1': '"#00BF86"',
                      'h3k27ac': 'purple',
                      'h3k27me3': 'chocolate',
                      'h4k16ac': 'purple',
                      'ha': 'orange',
                      'igg': 'grey',
                      'input': 'grey',
                      'normal': 'grey',
                      'oct4': 'gold',
                      'p53': 'gold',
                      'pbrm1': '"#9F2A90"',
                      'pou2f3': '#9F2A90',
                      'rad21': '"#910024"',
                      'rpb1': '#6FFF00',
                      's1184a': 'gold',
                      's1184e': 'red',
                      'smarca2': 'teal',
                      'smarca4': 'blue',
                      'smarcb1': '"#B8D1FF"',
                      'smarcc1': 'green',
                      'smarcd1': '"#BF49A6"',
                      'smarcd2': '"#BF49A6"',
                      'ss18': '"#6FFF00"',
                      'wt': 'blue'}


def get_colour(mark: str, mark_palette: Dict[str, str] = DEFAULT_COLOUR_MAP) \
        -> str:
    """
    Gets a standardised colour for the given mark.
    :param mark: The mark of interest.
    :param mark_palette: Palette to be used.
    """
    if mark not in mark_palette:
        raise ValueError(f'{mark} not found in palette.')

    return mark_palette[mark]
