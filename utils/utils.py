from copy import copy
from pathlib import Path
from typing import List, Dict, Union


def combine_cmds(cmds: List[str], num_per: int) -> List[str]:
    """
    Divides the cmds into (len(cmds) // num_per + 1) combined cmds,
    where each cmd contains at most <num_per> of the original commands.
    :param cmds: A list of commands.
    :param num_per: The number of cmds to assign to each job.
    :return: A list of jobs that will execute all cmds.
    """
    combined_cmds = []
    to_exec = []
    for cmd in cmds:
        to_exec.append(cmd)
        if len(to_exec) == num_per:
            combined_cmds.append('\n'.join(to_exec))
            to_exec = []
    if to_exec:
        combined_cmds.append('\n'.join(to_exec))
    return combined_cmds


class ExperimentalDesign:
    """
    === Description ===
    Stores useful information regarding the design of the experiment.

    === Private Attributes ===

    samples: All samples used in the experiment.
    sample_to_condition: Maps a given sample to the treatment/condition.
    condition_to_samples: Maps a condition to all samples that were generated
        under it.
    conditions: The conditions in the experiments.
    """

    _samples: List[str]
    _sample_to_condition: Dict[str, str]
    _sample_to_rep_number: Dict[str, int]
    _condition_to_samples: Dict[str, List[str]]
    _conditions: List[str]

    def __init__(self, sample_to_condition: Dict[str, str],
                 sample_to_rep_number: Dict[str, int] = None) -> None:
        """
        Initialises this experimental design using the given sample to
        condition mapping.
        :param sample_to_condition: Maps a given sample to the
            condition/treatment under which it was generated.
        """
        self._samples = sorted(sample_to_condition.keys())
        self._sample_to_condition = copy(sample_to_condition)
        self._condition_to_samples = {}

        # Invert sample_to_condition dict.
        for condition in self._sample_to_condition.values():
            self._condition_to_samples[condition] = []
            for sample in self._samples:
                if self._sample_to_condition[sample] == condition:
                    self._condition_to_samples[condition].append(sample)

        self._conditions = list(self._condition_to_samples.keys())[:]

        # Arbitrarily assign rep numbers if none are given.
        if not sample_to_rep_number:
            self._sample_to_rep_number = {}
            for cond in self._condition_to_samples:
                cts_sorted = sorted(self._condition_to_samples[cond])
                for i, sample in enumerate(cts_sorted):
                    self._sample_to_rep_number[sample] = i + 1
        else:
            self._sample_to_rep_number = copy(sample_to_rep_number)

    def get_condition(self, sample: str) -> str:
        return self._sample_to_condition[sample]

    def get_samples(self, condition: str = None) -> List[str]:
        if not condition:
            return copy(self._samples)
        else:
            return self._condition_to_samples[condition]

    def get_conditions(self) -> List[str]:
        return copy(self._conditions)

    def get_rep_num(self, sample: str) -> int:
        return self._sample_to_rep_number[sample]


def write_samples_file(filenames: List[Union[Path, str]],
                       out_path: Path) -> None:
    """
    Given some filenames, generates a samples file.

    >>> filenames = ['20230811_CGRNA087_HCC44_Rep1_S15_R1_001.fastq']
    >>> write_samples_file(filenames, Path('samples.txt'))
    >>> with open('samples.txt', 'r') as inf:
    >>>     print(inf.readline())
    CGRNA087	20230811_CGRNA087_HCC44_Rep1_S15_R1_001


    :param filenames: A list of filenames, with the sample at the second
        position after splitting by '_'.
    :param out_path: Path to the text file to write.
    """
    parsed_filenames = []
    for filename in filenames:
        if isinstance(filename, Path):
            parsed_filenames.append(filename.name)
        else:
            parsed_filenames.append(filename)
    filenames = parsed_filenames

    lines = {}
    for filename in filenames:
        sample_id = filename.split('_')[1]
        name = filename.split('.')[0]
        lines[name] = sample_id

    with open(out_path, 'w') as out_file:
        for filename in sorted(lines):
            out_file.write(filename + '\t' + lines[filename] + '\n')


