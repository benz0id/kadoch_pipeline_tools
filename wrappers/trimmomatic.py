from pathlib import Path
from typing import Dict
from constants.program_paths  import progs
from utils.path_manager import cmdify
from wrappers.wrapper import ProgramWrapper


class Trimmomatic(ProgramWrapper):
    def _get_dependencies(self) -> Dict[str, str]:
        return {"java": "jdk-1.8u112",
                "trimmomatic": "0.36"}

    def trim_paired_end(self, r1_path: Path, r2_path: Path, threads: int,
                        r1_out: Path = None, r2_out: Path = None,
                        paired_out: Path = None) -> None:
        cmd = cmdify(
            'java -jar', progs.trimmomatic,
            'PE',
            '-threads', threads,
            r1_path, r2_path, r1_out, paired_out, r2_out,

        )

