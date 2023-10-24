from pathlib import Path
from typing import List, Tuple

from utils.path_manager import cmdify
from utils.utils import TargetedDesign

DEFAULT_SUBSTRACT = '--operation subtract ' \
                    '-bs 50 ' \
                    '--smoothLength 600 ' \
                    '--samFlagInclude 64'


def subtract_background(bams: List[Path],
                        design: TargetedDesign,
                        out_dir: Path,
                        args: str = DEFAULT_SUBSTRACT,
                        control_id: str = 'IgG',
                        num_cores: int = 1) -> Tuple[List[str], List[Path]]:

    bam_to_control = design.get_file_to_control(bams, control_id)

    cmds = []
    outfiles = []

    for bam in bam_to_control:
        control = bam_to_control[bam]

        control_sample = control.name.split('_')[1]
        outfile = out_dir / (bam.name[:-4] + '.sub' + control_sample + '.bw')

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
