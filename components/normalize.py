from pathlib import Path
from typing import List

from utils.utils import TargetedDesign


def normalise_bams(bams: List[Path],
                   design: TargetedDesign,

                   control_id: str = 'IgG') -> List[Path]:
    bam_to_control = {}

    for bam in bams:
        sample = bam.name.split('_')[1]
        treatment = design.get_treatment(sample)
        control = design.find_in_files(bams,
                                       treatments=[treatment],
                                       marks=[control_id], num_expected=1)[0]

