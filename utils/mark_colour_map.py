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

DEFAULT_COLOUR_MAP = {
    'arid1a': 'orange',
    'ha': 'orange',
    'arid1b': '"#FF8C00"',
    'brg1': 'blue',
    'smarca4': 'blue',
    'smarca2': 'teal',
    'smarcc1': 'green',
    'c1': '"#0091FF"',
    'h3k27ac': 'purple',
    'h4k16ac': 'purple',
    'igg': 'grey',
    'input': 'grey',
    'empty': 'grey',
    'wt': 'blue',
    's1184a': 'gold',
    's1184e': 'red',
    '9ala': 'green',
    'p53': 'gold',
    'brd9': 'green',
    'brd4': 'green',
    'oct4': 'gold',
    'bicra': '"#00BF86"',
    'gltscr1': '"#00BF86"',
    'bicral': '"#00BF86"',
    'arid2': 'gold',
    'ss18': '"#6FFF00"',
    'smarcd1': '"#BF49A6"',
    'smarcd2': '"#BF49A6"',
    'rad21': '"#910024"',
    'smarcb1': '"#B8D1FF"',
    'ctcf': '"#BF49A6"',
    'atac': 'red',
    'normal': 'grey',
    'gltscr': 'green',
    'h3k27me3': 'chocolate',
    'pbrm1': '"#9F2A90"',
    'pou2f3': "#9F2A90",
    'rpb1': "#6FFF00"
}


def get_colour(mark: str, mark_pallete: Dict[str, str] = DEFAULT_COLOUR_MAP) \
        -> str:
    """
      Gets a standardised colour for the given mark.
      :param mark: The mark of interest.
      :param mark_pallete: Palette to be used.
      """
    if mark not in mark_pallete:
        raise ValueError(f'{mark} not found in palette.')

    return mark_pallete[mark]
