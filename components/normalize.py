from pathlib import Path
from typing import List, Tuple

from utils.path_manager import cmdify
from utils.utils import TargetedDesign

DEFAULT_SUBSTRACT = '--subtract -bs 50 --smoothLength 600 --samFlagInclude 64'

def subtract_background(bams: List[Path],
                        design: TargetedDesign,
                        out_dir: Path,
                        args: str = DEFAULT_SUBSTRACT,
                        control_id: str = 'IgG',
                        num_cores: int = 1) -> Tuple[List[str], List[Path]]:
    bam_to_control = {}

    for bam in bams:
        sample = bam.name.split('_')[1]
        treatment = design.get_treatment(sample)
        control = design.find_in_files(bams,
                                       treatments=[treatment],
                                       marks=[control_id], num_expected=1)[0]
        bam_to_control[bam] = control

    cmds = []
    outfiles = []

    for bam in bam_to_control:
        control = bam_to_control[bam]

        control_sample = control.name.split('_')[1]
        outfile = out_dir / (bam.name[:-4] + 'sub' + control_sample + '.bw')

        cmd = cmdify(
            "bamCompare",
            '--bamfile1', bam,
            '--bamfile2', control,
            '-o', outfile,
            '--numberOfProcessors', num_cores,
            args
        )

        cmds.append(cmd)
        outfiles.append(outfile)

    return cmds, outfiles





