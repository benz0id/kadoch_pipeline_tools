from pathlib import Path
from typing import List

from utils.utils import TargetedDesign


def normalise_bams(bams: List[Path],
                   design: TargetedDesign) -> List[Path]:

