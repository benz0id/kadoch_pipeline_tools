from copy import copy
from pathlib import Path
from typing import List, Dict, Union

from sample_sheet import SampleSheet


class DesignError(Exception):
    pass


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
    _sample_to_sample_id: Dict[str, str]
    _condition_to_samples: Dict[str, List[str]]
    _conditions: List[str]

    def __init__(self, sample_to_condition: Dict[str, str],
                 sample_to_rep_number: Dict[str, int] = None,
                 sample_to_sample_id: Dict[str, str] = None) -> None:
        """
        Initialises this experimental design using the given sample to
        condition mapping.
        :param sample_to_condition: Maps a given sample to the
            condition/treatment under which it was generated.
        """
        self._samples = sorted(sample_to_condition.keys())
        self._sample_to_condition = copy(sample_to_condition)
        self._condition_to_samples = {}
        self._sample_to_sample_id = sample_to_sample_id

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

    def get_invalid_sample_inds(self, strs: List[str],
                                fail_on_duplication: bool = True) -> List[int]:
        """
        Returns the indices of all elements of <strs> that do not contain a sample
        id as one of their substrings. Raises and error if a single sample is
        present more than once as a substring of some <strs>.
        :param fail_on_duplication:
        :param strs: A list of strings to be validated.
        :return: The indices of all invalid <strs>
        """

        samples = copy(self._samples)
        seen_samples = []
        invalid_samples = []

        for i, s in enumerate(strs):
            sample_found = False

            for sample in samples:
                if fail_on_duplication and \
                        sample in s and \
                        sample in seen_samples:
                    raise DesignError(sample + ' duplicated in ' + str(strs))

                if sample in s:
                    seen_samples.append(seen_samples)
                    sample_found = True

            if not sample_found:
                invalid_samples.append(i)

        return invalid_samples


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
        lines[sample_id] = name

    with open(out_path, 'w') as out_file:
        for filename in sorted(lines):
            out_file.write(filename + '\t' + lines[filename] + '\n')


def get_model_from_sample_sheet(sample_sheet: Path,
                                condition_slice: slice = slice(3, -1),
                                verbose: bool = False) -> ExperimentalDesign:
    """
    Infer the experimental model from the sample sheet.


    Expects sample names of the format:
    <date>_<sample_name>_<cell_type>_<conditions...>_Rep<rep_number>.

    :param sample_sheet: Path to the sample sheet.
    :param condition_slice: After splitting the sample id byu '_',
        the indices of the components that contain sample info.
    :param verbose: Whether to print out experimental design info.
    :return:
    """
    sample_sheet = SampleSheet(sample_sheet)
    sample_to_condition = {}
    sample_to_rep_number = {}
    sample_to_sample_id = {}

    for sample in sample_sheet.samples:
        sample_id = str(sample.sample_id)
        sample_name = sample_id.split('_')[1]

        if 'rep' in sample_id.lower():
            replicate_num = int(sample_name.lower().split('rep')[1])
        else:
            replicate_num = None

        condition = '_'.join(sample_id.split('_')[condition_slice])

        sample_to_condition[sample_name] = condition
        sample_to_rep_number[sample_name] = replicate_num
        sample_to_sample_id[sample_name] = sample_id

    if any([val is None
            for val in sample_to_rep_number.values()]):
        sample_to_rep_number = None

    return ExperimentalDesign(sample_to_condition,
                              sample_to_rep_number,
                              sample_to_sample_id)

















