from copy import copy
from pathlib import Path
from typing import List, Dict, Union, Any

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

    def get_reps(self) -> List[int]:
        return [self.get_rep_num(sample) for sample in self._samples]

    def get_sample_descs(self) -> List[str]:
        """
        :return:  a list of sample descriptors, formatted
            <condition>_Rep<rep_num>.
        For sample in samples. Is in the same order as samples.
        """
        descs = []

        # If there are no replicates, just return to which group each sample
        # belongs.
        if len(set(self._sample_to_rep_number.values())) == 1:
            for sample in self._samples:
                descs.append(self._sample_to_condition[sample])
            return descs
        # Append replicate number to end of description.
        else:
            for sample in self._samples:
                cond = self._sample_to_condition[sample]
                rep = self.get_rep_num(sample)
                descs.append(cond + '_Rep' + str(rep))
            return descs

    def align_files_to_samples(self, files: List[Path]) -> List[Path]:
        """
        Given that each filename contains a sample name, allowing for a 1 to 1
        mapping between each of files and samples, returns files ordered by
        sample.
        :param files: A list of filepaths each containing sample names.
        :return: The same list, ordered in the same way as samples.
        """
        rtrn = []

        for sample in self._samples:
            match_found = False

            for file in files:
                if sample in file.name and match_found:
                    raise ValueError(f"Multiple matches found for {sample}")
                elif sample in file.name:
                    match_found = True
                    rtrn.append(file)

            if not match_found:
                raise ValueError(f"Could not find match for {sample}")

        return rtrn

    def merge_reps(self):
        """
        Returns the design of the experiment after merging replicates.
        :return:
        """

        sample_to_condition = {}
        sample_to_rep_number = {}
        sample_to_sample_id = {}

        for condition in self._conditions:
            samples = sorted(self.get_samples(condition))
            merged_sample_name = '-'.join(samples)

            date_string = self._sample_to_sample_id[samples[0]][:8]
            merged_id = f'{date_string}_{merged_sample_name}_{condition}'
            sample_to_condition[merged_sample_name] = condition
            sample_to_rep_number[merged_sample_name] = 0
            sample_to_sample_id[merged_sample_name] = merged_id

        return ExperimentalDesign(sample_to_condition,
                                  sample_to_rep_number,
                                  sample_to_sample_id)

    def get_rep_num(self, sample: str) -> int:
        return self._sample_to_rep_number[sample]

    def __str__(self) -> str:
        def pretty_dict(dic: Dict[str, Any]) -> str:
            s = ''
            for sample in self._samples:
                s += '\t' + sample + ': ' + str(dic[sample]) + '\n'
            return s

        return \
            f'=== Sample: Sample ID ===\n' \
            f'{pretty_dict(self._sample_to_sample_id)}\n' \
            f'=== Sample: Condition ===\n' \
            f'{pretty_dict(self._sample_to_condition)}\n' \
            f'=== Sample: Replicate ===\n' \
            f'{pretty_dict(self._sample_to_rep_number)}\n'

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

    def get_fastq_groupings(self, fastqs: List[Path]) -> List[List[Path]]:
        """
        Groups the fastqs according to which merged sample they belong to.
        :return:
        """
        pe_found = False
        se_found = False
        for file in fastqs:
            tail = file.name.split('.')[0].split('_')[-2]
            if "R1" in tail or "R2" in tail:
                pe_found = True
            else:
                se_found = True

        if pe_found and se_found:
            raise ValueError("A mixture of paired end and single end reads"
                             "were found.")

        if pe_found:
            r1s = []
            r2s = []
            for file in fastqs:
                tail = file.name.split('.')[0].split('_')[-2]
                if "R1" in tail:
                    r1s.append(file)
                elif "R2" in tail:
                    r2s.append(file)
                else:
                    raise ValueError(
                        "I have no idea how this is even remotely "
                        "possible. Sorry :)")

            r1_groups = self.get_groupings(r1s)
            r2_groups = self.get_groupings(r2s)

            return r1_groups + r2_groups
        else:
            return self.get_groupings(fastqs)

    def get_groupings(self, files: List[Path]) -> List[List[Path]]:
        """
        For experimental designs that have been merged with others. Returns the
        groups of samples that have been merged.

        Useful for merging fastqs and bam files.

        :param files: A list of files. For each sample contained in this
            design, <files> must contain exactly one file with the sample's
            name in it.
        :return: Files grouped by the merged sample they belong to.
        """

        sample_groups = [sample.split('-') for sample in self._samples]
        groups = [[] for _ in range(len(sample_groups))]

        for file in files:
            found = False
            for i, sample_group in enumerate(sample_groups):
                for sample in sample_group:
                    match = sample in file.name
                    if match and found:
                        raise ValueError("Multiple matches found for" +
                                         str(file))
                    elif match:
                        found = True
                        groups[i].append(file)
            if not found:
                raise ValueError("Could not find match for" + str(file))

        return groups


def combine_runs(*runs: ExperimentalDesign) -> ExperimentalDesign:
    """
    Combines the give experimental designs into a single experimental design.
    All samples with matching condition and replicate number are merged into
    a single sample.

    :param runs: A list of experiments.
    :return: The merged experiment.
    """

    cond_rep_to_samples = {}

    def add(key, val):
        if key in cond_rep_to_samples:
            cond_rep_to_samples[key].append(val)
        else:
            cond_rep_to_samples[key] = [val]

    for run in runs:
        for sample in run.get_samples():
            cond = run.get_condition(sample)
            rep = run.get_rep_num(sample)
            add((cond, rep), sample)

    sample_to_condition = {}
    sample_to_rep_number = {}
    sample_to_sample_id = {}

    for cond, rep in cond_rep_to_samples:
        samples = cond_rep_to_samples[(cond, rep)]
        sample = '-'.join(samples)
        sample_id = '_'.join(['0MERGED0', sample, cond, 'Rep' + str(rep)])

        sample_to_condition[sample] = cond
        sample_to_rep_number[sample] = rep
        sample_to_sample_id[sample] = sample_id

    return ExperimentalDesign(sample_to_condition,
                              sample_to_rep_number,
                              sample_to_sample_id)


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
                                condition_slice: slice = slice(2, -1),
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
            replicate_num = int(sample_id.lower().split('rep')[1])
        else:
            replicate_num = None

        condition = '_'.join(sample_id.split('_')[condition_slice])

        sample_to_condition[sample_name] = condition
        sample_to_rep_number[sample_name] = replicate_num
        sample_to_sample_id[sample_name] = sample_id

    if any([val is None
            for val in sample_to_rep_number.values()]):
        sample_to_rep_number = None

    if verbose:
        print(len(sample_to_condition), sample_to_condition, '\n')
        print(len(sample_to_rep_number), sample_to_rep_number, '\n')
        print(len(sample_to_sample_id), sample_to_sample_id, '\n')

    return ExperimentalDesign(sample_to_condition,
                              sample_to_rep_number,
                              sample_to_sample_id)
