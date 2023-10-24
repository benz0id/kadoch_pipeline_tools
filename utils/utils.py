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

    def remove_str(self, s: str, samples: List[str] = None) -> None:
        """
        Removes the given string from all attributes associated with the
        given samples.
        :param s: The string to be removed
        :param samples: The samples to remove the given string from. All
            samples by default.
        """
        if not samples:
            samples = self._samples

        conditions = []
        for sample in samples:
            conditions.append(self._sample_to_condition[sample])

        while conditions:
            condition = conditions.pop()
            pruned_condition = condition.replace(s, '')
            pruned_condition = pruned_condition.replace('__', '_')

            if pruned_condition[0] == '_':
                pruned_condition = pruned_condition[1:]

            if pruned_condition[-1] == '_':
                pruned_condition = pruned_condition[:-1]

            if pruned_condition in self._conditions:
                continue
            else:
                print(condition, pruned_condition)
                for sample in self._condition_to_samples[condition]:
                    self._sample_to_condition[sample] = pruned_condition

                ind = self._conditions.index(condition)
                self._conditions.remove(condition)
                self._conditions.insert(ind, pruned_condition)

                self._condition_to_samples[pruned_condition] = \
                    self._condition_to_samples[condition]

                self._condition_to_samples.pop(condition)

    def get_condition(self, sample: str) -> str:
        return self._sample_to_condition[sample]

    def get_samples(self, condition: str = None,
                    rep: int = None) -> List[str]:

        matching_samples = []
        for sample in self._samples:

            if condition and self._sample_to_condition[sample] == condition:
                continue
            if rep and not self._sample_to_rep_number[sample] == rep:
                continue
            matching_samples.append(sample)

        return matching_samples

    def get_conditions(self) -> List[str]:
        return copy(self._conditions)

    def get_sample_id(self, sample: str) -> str:
        return self._sample_to_sample_id[sample]

    def get_ordered_conditions(self) -> List[str]:
        return [self.get_condition(sample) for sample in self._samples]

    def get_ordered_reps(self) -> List[int]:
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

    def align_to_samples(self, to_align: Union[List[str], List[Path]]) \
            -> Union[List[Path], List[str]]:
        """
        Given that each sample_str contains a sample name, allowing for a 1 to 1
        mapping between each of sample_strs and samples, returns files ordered by
        sample.
        :param to_align: A list of sample_strs each containing sample names.
        :return: The same list, ordered in the same way as samples.
        """
        return align_to_strs(to_align, self._samples)

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
            date_string = self._sample_to_sample_id[samples[0]][:8]

            # Split sample that have already been merged to maintain order.
            s_split = []
            for sample in samples:
                s_split.extend(sample.split('-'))
            samples = s_split

            samples.sort()
            merged_sample_name = '-'.join(samples)
            merged_id = f'{date_string}_{merged_sample_name}_{condition}'
            sample_to_condition[merged_sample_name] = condition
            sample_to_rep_number[merged_sample_name] = 0
            sample_to_sample_id[merged_sample_name] = merged_id

        return ExperimentalDesign(sample_to_condition,
                                  sample_to_rep_number,
                                  sample_to_sample_id)

    def get_rep_num(self, sample: str) -> int:
        return self._sample_to_rep_number[sample]

    def remove_condition(self, conditions: Union[str, List[str]]) -> None:
        """
        Remove all sample corresponding to the given conditions. All
        <conditions> must be contained within this design.
        :return: None
        """

        if isinstance(conditions, str):
            conditions = [conditions]

        for condition in conditions:
            samples = self.get_samples(condition)
            for sample in samples:
                self.remove_sample(sample)

            del self._condition_to_samples[condition]
            self._conditions.remove(condition)

    def remove_sample(self, sample: str) -> None:
        self._samples.remove(sample)
        del self._sample_to_condition[sample]
        del self._sample_to_rep_number[sample]
        del self._sample_to_sample_id[sample]

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

    def get_fastq_groupings(self, fastqs: List[Path], verbose=False) -> List[
        List[Path]]:
        """
        Groups the fastqs according to which merged sample they belong to.

        Each group at return[i] corresponds to the merged samples at
        self.get_samples[i]

        If paired end inputs are received, the returned list will

        :return:
        """
        f_found = False
        r_found = False
        for file in fastqs:
            tail = file.name.split('.')[0].split('_')[-2]
            if "R1" in tail:
                f_found = True
            elif "R2" in tail:
                r_found = True

        if f_found and r_found:
            if verbose:
                print('Paired end mode')

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

            r1_groups = self.get_groupings(r1s, verbose=verbose)
            r2_groups = self.get_groupings(r2s, verbose=verbose)

            if verbose:
                for i, group in enumerate(r1_groups):
                    s = str([f.name.split('_')[1] for f in group])
                    print(s, '-',
                          self.get_ordered_conditions()[i], '-',
                          self.get_ordered_reps()[i])

            return r1_groups + r2_groups

        else:
            return self.get_groupings(fastqs)

    def get_groupings(self, files: List[Path], verbose: bool = False) \
            -> List[List[Path]]:
        """
        For experimental designs that have been merged with others. Returns the
        groups of samples that have been merged.

        Useful for merging fastqs and bam files.

        :param verbose:
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
                added = False

                for sample in sample_group:
                    match = sample in file.name

                    if verbose and match and not found:
                        print(file.name, '->', sample_group)

                    if match and found:
                        raise ValueError("Multiple matches found for" +
                                         str(file))
                    elif match and not added:
                        groups[i].append(file)
                        added = True

                if added:
                    found = True
            if verbose:
                print('\n')
            if not found:
                raise ValueError("Could not find match for" + str(file))

        return groups

    def find_in_files(self, files: List[Path],
                      samples: List[str] = None,
                      conditions: List[str] = None,
                      reps: List[int] = None) -> List[Path]:
        valid_files = []
        for file in files:
            sample = file.name.split('_')[1]

            if sample not in self._samples:
                raise ValueError(
                    f'{sample} in {file.name} is not a valid sample.')

            cond = self._sample_to_condition[sample]
            rep = self._sample_to_rep_number[sample]

            valid_sample = not samples or sample in samples
            valid_cond = not conditions or cond in conditions
            valid_rep = not reps or rep in reps

            if valid_rep and valid_cond and valid_sample:
                valid_files.append(file)
        return valid_files


class TargetedDesign(ExperimentalDesign):

    def __init__(self, sample_to_condition: Dict[str, str],
                 sample_to_mark: Dict[str, str],
                 sample_to_treatment: Dict[str, str],
                 sample_to_rep_number: Dict[str, int] = None,
                 sample_to_sample_id: Dict[str, str] = None) -> None:
        super().__init__(sample_to_condition, sample_to_rep_number,
                         sample_to_sample_id)

        self._sample_to_mark = copy(sample_to_mark)
        self._sample_to_treatment = copy(sample_to_treatment)

    def get_samples(self, condition: str = None,
                    mark: str = None,
                    treatment: str = None,
                    rep: int = None) -> List[str]:

        matching_samples = []
        for sample in super().get_samples(condition, rep):
            if mark and not self._sample_to_mark[sample] == mark:
                continue
            if treatment and not self._sample_to_treatment[
                                     sample] == treatment:
                continue
            matching_samples.append(sample)

        return matching_samples

    def get_mark(self, sample: str) -> str:
        return self._sample_to_mark[sample]

    def get_file_to_control(self, files: List[Path],
                            control_id: str = 'IgG') -> Dict[Path, Path]:
        file_to_control = {}
        for bam in files:
            sample = bam.name.split('_')[1]
            treatment = self.get_treatment(sample)
            control = self.find_in_files(files,
                                         treatments=[treatment],
                                         marks=[control_id], num_expected=1)[
                0]
            file_to_control[bam] = control
        return file_to_control

    def get_marks(self, samples: List[str] = None) -> List[str]:
        if not samples:
            samples = self.get_samples()
        return [self.get_mark(sample) for sample in samples]

    def get_treatment(self, sample: str) -> str:
        return self._sample_to_treatment[sample]

    def get_treatments(self, samples: List[str] = None) -> List[str]:
        if not samples:
            samples = self.get_samples()
        return [self.get_treatment(sample) for sample in samples]

    def find_in_files(self, files: List[Path],
                      samples: List[str] = None,
                      conditions: List[str] = None,
                      reps: List[int] = None,
                      marks: List[str] = None,
                      treatments: List[str] = None,
                      num_expected: int = None) -> List[Path]:
        valid_files = []

        super_valid = super().find_in_files(files=files, samples=samples,
                                            conditions=conditions, reps=reps)
        for file in super_valid:
            sample = file.name.split('_')[1]

            if sample not in self._samples:
                raise ValueError(
                    f'{sample} in {file.name} is not a valid sample.')
            mark = self.get_mark(sample)
            treat = self.get_treatment(sample)

            valid_mark = not marks or mark in marks
            valid_treat = not treatments or treat in treatments

            if valid_mark and valid_treat:
                valid_files.append(file)

        if num_expected and len(valid_files) != num_expected:
            raise ValueError("More files found than specified")

        return valid_files

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
            f'{pretty_dict(self._sample_to_rep_number)}\n' \
            f'=== Sample: Mark ===\n' \
            f'{pretty_dict(self._sample_to_mark)}\n' \
            f'=== Sample: Treatment ===\n' \
            f'{pretty_dict(self._sample_to_treatment)}\n'

    def remove_sample(self, sample: str) -> None:
        super().remove_sample(sample)
        del self._sample_to_mark[sample]
        del self._sample_to_treatment[sample]


def convert_to_targeted(design: ExperimentalDesign, mark_slice: slice,
                        treatment_slice: slice) -> TargetedDesign:
    """
    Converts the given experimental design into a targeted design.
    :param design: An experimental design. Cannot be an instance of Targetted
        design.
    :return: TargetedDesign
    """
    samples = design.get_samples()

    sample_to_cond = {}
    cond_to_samples = {}
    sample_to_rep_num = {}
    sample_to_mark = {}
    sample_to_treat = {}
    sample_to_sample_id = {}

    def add(dic, key, val) -> None:
        if key in dic:
            dic[key].append(val)
        else:
            dic[key] = [val]

    for sample in samples:
        cond = design.get_condition(sample)
        sample_to_cond[sample] = cond
        add(cond_to_samples, cond, sample)
        sample_to_rep_num[sample] = design.get_rep_num(sample)
        mark = '_'.join(cond.split('_')[mark_slice])
        treat = '_'.join(cond.split('_')[treatment_slice])
        print(cond.split('_'), mark, treat)
        sample_to_mark[sample] = mark
        sample_to_treat[sample] = treat

        sample_to_sample_id[sample] = design.get_sample_id(sample)

    return TargetedDesign(sample_to_cond,
                          sample_to_mark,
                          sample_to_treat,
                          sample_to_rep_num,
                          sample_to_sample_id)


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

        # Split sample that have already been merged to maintain order.
        s_split = []
        for sample in samples:
            s_split.extend(sample.split('-'))
        samples = s_split

        samples.sort()
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


def align_to_strs(to_align: Union[List[str], List[Path]],
                  strs: List[str]) \
        -> Union[List[Path], List[str]]:
    """
    Given that each sample_str contains a sample name, allowing for a 1 to 1
    mapping between each of to_align and strs, returns files ordered by
    strs.

    Note: There must be a one to one mapping between <strs> and <to_align>,
        where each of <to_align> contains one of <strs> as a substring.

    If to_align is a list of filepaths, the filenames will be considered for
    ordering.

    :param to_align: A list of strings or paths.
    :param strs: A list of strings.
    :return: The same list, ordered in the same way as samples.
    """
    rtrn = []

    sample_to_matches = {}
    matches_to_sample = {}
    order = []

    def add(d: Dict, key: Any, val: Any):
        if key in d:
            d[key].append(val)
        else:
            d[key] = [val]

    # Find order and mappings.
    for sample in strs:
        for i, p_to_align in enumerate(to_align):
            # Convert from path to string if necessary.
            if isinstance(p_to_align, Path):
                s_to_align = p_to_align.name
            elif isinstance(p_to_align, str):
                s_to_align = p_to_align
            else:
                ValueError("Unknown Input Type")

            if sample in s_to_align:
                order.append(i)
                add(sample_to_matches, sample, p_to_align)
                add(matches_to_sample, p_to_align, sample)

    # Ensure one to one mappings.
    for sample in strs:
        if sample not in sample_to_matches:
            raise ValueError(f"Failed to find match for {sample}.")

        if len(sample_to_matches[sample]) > 1:
            raise ValueError(f"Found multiple matches for {sample}: "
                             f"{str(sample_to_matches[sample])}.")

    for to_match in to_align:
        if to_match not in matches_to_sample:
            raise ValueError(f"Failed to find sample for {to_match}.")
        if len(matches_to_sample[to_match]) > 1:
            raise ValueError(f"Found multiple sample matches for "
                             f"{to_match}: {str(matches_to_sample[to_match])}.")

    ordered = []
    for ind in order:
        ordered.append(to_align[ind])

    return ordered
