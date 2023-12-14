import dataclasses
import sys
from copy import copy, deepcopy
from pathlib import Path
from typing import List, Dict, Union, Any, Tuple

from sample_sheet import SampleSheet


class DesignError(Exception):
    pass


def check_one_to_one(a: List[str], b: List[str], name: str) -> None:
    """
    Checks whether a and b have a one to one mapping.
    :param a: A list of strings.
    :param b: A list of strings
    :param name: Common feature between the mappings.
    :return:
    """
    unrepresented = copy(a)
    # All conditions in this sample must be represented in this condition.
    for s in b:
        if s not in a:
            raise ValueError(f'{s} not in the provided '
                             f'{name} order list.')
        if s in unrepresented:
            unrepresented.remove(s)

    if unrepresented:
        unrepr_str = ', '.join(unrepresented)
        print(f'The following conditions were not found in this '
              f'design: {unrepr_str}', file=sys.stderr)

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


MANDATORY_SAMPLE_FIELDS = ['sample_name', 'condition', 'replicate']
TARGETED_SAMPLE_FIELDS = ['target', 'treatment']
KNOWN_OPTIONAL_FIELDS = ['bam', 'bai', 'bed', 'narrow_peak', 'broad_peak',
                         'bw', 'big_wig']


@dataclasses.dataclass
class Sample:
    """
    === Description ===
    Stores attributes found in all samples.

    === Public Attributes ===
    sample_name: The ID of the sample. Usually the experimenters initials
        followed by the protocol and then the sample number.
    condition: All information regarding the conditions used to generate this
        sample.
    replicate: The replicate number of this sample.
    """
    sample_name: str
    condition: str
    replicate: int

    _condition_order: List[str]

    target: str
    treatment: str

    bam: Path
    bai: Path
    bed: Path
    narrow_peak: Path
    broad_peak: Path

    def __init__(self, sample_name: str, condition: str,
                 replicate: int, **kwargs) -> None:
        """
        :param sample_name: The ID of the sample. Usually the experimenters
            initials followed by the protocol and then the sample number.
        :param condition: All information regarding the conditions used to
            generate this sample.
        :param replicate: The replicate number of this sample.
        """
        self.sample_name = sample_name
        self.condition = condition
        self.replicate = replicate

        for arg in kwargs:
            self.__dict__[arg] = kwargs[arg]

    def get_attrs(self) -> Dict[str, Any]:
        """
        Return all of this sample's attributes.
        :return:
        """
        return self.__dict__

    def update_attrs(self, attrs: Dict[str, Any]) -> None:

        for attr in attrs:
            val = attrs[attr]

            if val is None:
                continue

            if attr not in self.__dict__:
                raise ValueError(f'Unrecognised attribute: {attr}')

    def __eq__(self, other) -> bool:
        return self.sample_name == other.sample_name

    def set_condition_order(self, condition_order: List[str]) -> None:
        """
        Sets the order of the conditions for internal ordering comparisons.
        :param condition_order: 
        :return: 
        """
        self._condition_order = condition_order

    def __lt__(self, other) -> bool:
        if '_condition_order' not in self.__dict__:
            return self.sample_name < other.sample_name

        assert self._condition_order == other._condition_order
        return self._condition_order.index(self.condition) < \
            self._condition_order.index(other.condition)


class TargetedSample(Sample):
    """
    === Description ===
    Stores attributes of targeted samples which target a certain mark, e.g.
    CUT&RUN, ChIP.

    === Public Attributes ===
    sample_name: The ID of the sample. Usually the experimenters initials
        followed by the protocol and then the sample number.
    condition: All information regarding the conditions used to generate this
        sample.
    replicate: The replicate number of this sample.
    target: The target of this sample.
    treatment: The treatment applied to this Sample. The same as
            <condition> in a generic sample.


    === Private Attributes ===
    precedence: Whether to sort first by target or treatment.
    target_order: A list of targets, in the order into which they should be
        sorted.
    treatment_order: A list of treatments, in the order into which they should
        be sorted.
    """
    target: str
    treatment: str

    _precedence: str
    _target_order: List[str]
    _treatment_order: List[str]

    def __init__(self, sample_name: str, condition: str, replicate: int,
                 target: str, treatment: str, **kwargs) -> None:
        """
        :param sample_name: The name of this sample.
        :param condition: The condition under which this sample was generated.
            Includes both the target/mark and the treatement applied.
        :param replicate: The replicate of this sample.
        :param target: The target of this sample.
        :param treatment: The treatment applied to this Sample. The same as
            <condition> in a generic sample.
        :param kwargs: All other attributes of this sample.
        """
        super().__init__(sample_name, condition, replicate, **kwargs)
        self.target = target
        self.treatment = treatment
        self._condition_order = []
        self._precedence = 'treatment'

    def set_target_order(self, target_order: List[str]) -> None:
        """
        Sets the order of the conditions for internal ordering comparisons.
        :param target_order:
        :return: 
        """
        self._target_order = target_order

    def set_treatment_order(self, treatment_order: List[str]) -> None:
        """
        Sets the order of the conditions for internal ordering comparisons.
        :param treatment_order: 
        :return: 
        """
        self._treatment_order = treatment_order

    def __lt__(self, other) -> bool:
        if '_target_order' not in self.__dict__ or \
                '_treatment_order' not in self.__dict__:
            return self.sample_name < other.sample_name

        assert self._target_order == other._target_order
        assert self._treatment_order == other._treatment_order

        self_targ = self._target_order.index(self.target)
        self_treat = self._treatment_order.index(self.treatment)

        other_targ = self._target_order.index(other.target)
        other_treat = self._treatment_order.index(other.treatment)

        if self._precedence == 'target':
            self_targ = self_targ * len(self._treatment_order)
            other_targ = other_targ * len(self._treatment_order)
        elif self._precedence == 'treatment':
            self_treat = self_treat * len(self._target_order)
            other_treat = other_treat * len(self._target_order)
        else:
            raise ValueError(f"{self._precedence} is not a valid precedence.")

        return self_treat + self_targ < other_treat + other_targ


def args_to_filters(**kwargs) -> Dict[str, Any]:
    """
    Converts a dictionary of arguments into a dictionary or filters.
    Used to adapt old implementation to new implementation.
    :param kwargs:
    :return:
    """

    # Add more as needed.
    valid_filtering_arguments = [
        'sample_name',
        'condition',
        'replicate',
        'target',
        'treatment'
    ]
    kas = ', '.join(valid_filtering_arguments)

    # Convert to filters.
    filters = {}
    for arg in kwargs:
        if kwargs[arg] is None:
            continue

        val = kwargs[arg]

        if arg == 'sample':
            arg = 'sample_name'

        if arg == 'marks' or arg == 'mark':
            arg = 'target'

        # Old convention was plural - new convention is singular
        if arg[-1] == 's':
            arg = arg[:-1]

        if arg not in valid_filtering_arguments:
            raise ValueError(f'Adapter function received unsupported '
                             f'sample attribute {arg}, known attributes '
                             f'include {kas}.')

        filters[arg] = val

    return filters


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
    _samples: List[Sample]
    _condition_order: List[str]

    def __init__(self, sample_to_condition: Dict[str, str],
                 sample_to_rep_number: Dict[str, int] = None) -> None:
        """
        Initialises this experimental design using the given sample to
        condition mapping.
        :param sample_to_condition: Maps a given sample to the
            condition/treatment under which it was generated.
        """
        cond_counts = {}
        self._samples = []

        for sample in sample_to_condition:
            condition = sample_to_condition[sample]

            # Manually assign replicate numbers if replicates aren't provided.
            if sample_to_rep_number:
                rep_num = sample_to_rep_number[sample]
            else:
                if condition in cond_counts:
                    cond_counts[condition] += 1
                    rep_num = cond_counts[condition]
                else:
                    cond_counts[condition] = 1
                    rep_num = cond_counts[condition]

            sample = Sample(sample, condition, rep_num)
            self._samples.append(sample)

    def set_condition_order(self, condition_order: List[str]) -> None:
        """
        Sets the order in which to sort the samples in this design.
        :param condition_order: The order of the sample conditions. Must 
            contain conditions for all samples.
        :return:
        """

        check_one_to_one(condition_order,
                         [sample.condition for sample in self._samples],
                         'condition')

        self._condition_order = condition_order

        for sample in self._samples:
            sample.set_condition_order(condition_order)

    def sort_samples(self) -> None:
        self._samples.sort()

    def remove_str(self, s: str, samples: List[str] = None, pos: int = None,
                   attrs: List[str] = None) -> None:
        """
        Removes the given string from all attributes associated with the
        given samples.
        :param s: The string to be removed
        :param samples: The samples to remove the given string from. All
            samples by default.
        :param pos: Only remove <s> if it occurs at <pos> after splitting by
        '_'.
        :param attrs: A list of all the attributes to modify. All attributes by
        default.
        """

        for sample in self._samples:

            if samples and sample.sample_name not in samples:
                continue

            sample_attrs = sample.get_attrs()

            for attr in sample_attrs:
                val = sample_attrs[attr]

                if attrs and attr not in attrs:
                    continue

                if not isinstance(val, str):
                    continue

                if pos is not None:
                    components = val.split('_')
                    if s == components[pos]:
                        components.pop(pos)
                    pruned_val = '_'.join(components)
                else:
                    pruned_val = val.replace(s, '')
                    pruned_val = pruned_val.replace('__', '_')

                if pruned_val[0] == '_':
                    pruned_val = pruned_val[1:]

                if pruned_val[-1] == '_':
                    pruned_val = pruned_val[:-1]

                sample_attrs[attr] = pruned_val
            sample.update_attrs(sample_attrs)

    def query(self, get: str,
              filters: Dict[str, Any] = None) -> List[Any]:
        """
        Returns <rtrn_attr> from each sample matching <filters>.
        :param filters: A dictionary mapping attribute names to one or more
            allowable values as a list.
        :param get: The attribute to return from all samples matching
            <filters>.
        :return: A list of attributes of type <rtrn_attrs>.
        """

        if filters is None:
            filters = {}

        filters = copy(filters)
        for filter in filters:
            if not isinstance(filters[filter], list):
                filters[filter] = [filters[filter]]

        valid_samples = []

        for sample in self._samples:
            attrs = sample.get_attrs()
            valid = True

            for attr in filters:

                if attr not in attrs:
                    allowable = ', '.join(list(attrs.keys()))
                    raise ValueError(
                        f'Received unexpected attribute {attr}. '
                        f'Allowble attributes include {allowable}.')

                allowable_vals = filters[attr]
                sample_val = attrs[attr]
                if sample_val not in allowable_vals:
                    valid = False
                    break

            if valid:
                valid_samples.append(sample)

        return_attrs = []
        for sample in valid_samples:
            attrs = sample.get_attrs()

            if get not in attrs:
                allowable = ', '.join(list(attrs.keys()))
                raise ValueError(
                    f'Received unexpected return attribute {get}. '
                    f'Allowble attributes include {allowable}.')

            return_attrs.append(attrs[get])
        return return_attrs

    def find_in_files(self, files: List[Path],
                      samples: List[str] = None,
                      conditions: List[str] = None,
                      reps: List[int] = None,
                      marks: List[str] = None,
                      treatments: List[str] = None,
                      num_expected: int = None) -> List[Path]:

        filters = args_to_filters(
            samples=samples,
            conditions=conditions,
            reps=reps,
            marks=marks,
            treatments=treatments
        )

        valid_samples = self.query('sample_name', filters)
        valid_files = []

        for file in files:
            sample_name = file.name.split('_')[1]

            if sample_name not in [sample.sample_name for sample in
                                   self._samples]:
                raise ValueError(
                    f'{sample_name} in {file.name} is not a valid sample.')

            if sample_name in valid_samples:
                valid_files.append(file)

        if num_expected and len(valid_files) != num_expected:
            f_str = '\t\n' + '\t\n'.join([f.name for f in valid_files])
            raise ValueError(f"Unexpected number of files found,"
                             f" {len(valid_files)} != {num_expected} {f_str}")

        return valid_files

    def get_samples(self, conditions: List[str] = None,
                    reps: List[int] = None,
                    num_expected: int = None) -> List[str]:
        filters = {}

        if conditions:
            filters['condition'] = conditions
        if reps:
            filters['reps'] = reps

        return self.query('sample_name', filters)

    def get_conditions(self) -> List[str]:
        return list(set([sample.condition for sample in self._samples]))

    def get_condition(self, sample_name: str) -> str:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                return sample.condition

    def get_ordered_conditions(self) -> List[str]:
        return [sample.condition for sample in self._samples]

    def get_ordered_reps(self) -> List[int]:
        return [self.get_rep_num(sample.sample_name) for sample in
                self._samples]

    def get_sample_descs(self) -> List[str]:
        """
        :return:  a list of sample descriptors, formatted
            <condition>_Rep<rep_num>.
        For sample in samples. Is in the same order as samples.
        """
        descs = []

        # If there are no replicates, just return to which group each sample
        # belongs.
        if len(set([sample.replicate for sample in self._samples])) == 1:
            for sample in self._samples:
                descs.append(sample.condition)
            return descs
        # Append replicate number to end of description.
        else:
            for sample in self._samples:
                cond = sample.condition
                rep = sample.replicate
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
        return align_to_strs(to_align,
                             [sample.sample_name for sample in self._samples])

    def merge_reps(self):
        """
        Returns the design of the experiment after merging replicates.
        :return:
        """

        sample_to_condition = {}
        sample_to_rep_number = {}

        for condition in self.get_conditions():
            samples = sorted(self.get_samples(conditions=[condition]))

            # Split sample that have already been merged to maintain order.
            s_split = []
            for sample in samples:
                s_split.extend(sample.split('-'))
            samples = s_split

            samples.sort()
            merged_sample_name = '-'.join(samples)
            sample_to_condition[merged_sample_name] = condition
            sample_to_rep_number[merged_sample_name] = 0

        return ExperimentalDesign(sample_to_condition,
                                  sample_to_rep_number)

    def get_rep_num(self, targ_sample_name: str) -> int:
        for sample in self._samples:
            if sample.sample_name == targ_sample_name:
                return sample.replicate

    def remove_condition(self, conditions: Union[str, List[str]]) -> None:
        """
        Remove all sample corresponding to the given conditions. All
        <conditions> must be contained within this design.
        :return: None
        """

        if isinstance(conditions, str):
            conditions = [conditions]
        to_rem = []
        for sample in self._samples:
            if sample.condition in conditions:
                to_rem.append(sample)

        for sample in to_rem:
            self._samples.remove(sample)

    def remove_sample(self, sample_name: str) -> None:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                self._samples.remove(sample)
                return
        raise ValueError(f'{sample_name} not found.')

    def __str__(self) -> str:
        def pretty_dict(dic: Dict[str, Any]) -> str:
            s = ''
            for sample in self._samples:
                s += '\t' + sample.sample_name + ': ' + \
                     str(dic[sample.sample_name]) + '\n'
            return s

        return \
            f'=== Sample: Condition ===\n' \
            f'{pretty_dict({sample.sample_name: sample.condition for sample in self._samples})}\n' \
            f'=== Sample: Replicate ===\n' \
            f'{pretty_dict({sample.sample_name: sample.replicate for sample in self._samples})}\n'

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

        samples = [sample.sample_name for sample in self._samples]
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
            tail = file.name.split('.')[0].split('_')[-2:]
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

        sample_groups = [sample.sample_name.split('-')
                         for sample in self._samples]
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

    def subset_design(self, matching: List[Tuple[int, str]] = None,
                      not_matching: List[Tuple[int, str]] = None,
                      include_default: bool = False):
        """
        Returns a subset of this design with the containing or not containing,
        the given strings at the positions given in <matching> and
        <not_matching>, respectively.
        :param matching: After splitting each condition by '_', must contain
            at least one of <str> at position <int> in order to be included.
        :param not_matching: After splitting each condition by '_', must not
            contain any of <str> at position <int> in order to be included.
        :param include_default: Whether samples should be included by default.
        :return: A design with the selected samples.
        """

        if not matching:
            matching = []

        if not not_matching:
            not_matching = []

        to_inc = []

        for sample in self._samples:
            components = sample.condition.split('_')
            inc = include_default

            for i, s in matching:
                if s in components[i]:
                    inc = True

            for i, s in not_matching:
                if s in components[i]:
                    inc = False

            if inc:
                to_inc.append(sample.sample_name)

        to_rem = set([sample.sample_name
                      for sample in self._samples]) - set(to_inc)
        subset = copy(self)
        for sample_name in to_rem:
            subset.remove_sample(sample_name)
        return subset

    def __add__(self, other):
        """
        Returns an experimental design containing samples from both
        experimental designs. Resets condition orders if they are not
        identical.
        :param other: Another experimental design.
        :return:
        """

        if type(self) != type(other):
            raise ValueError('Addition not defined between these classes.')

        rtrn = deepcopy(self)

        rtrn._samples.extend(other._samples)

        for cond in other._condition_order:
            if cond not in self._condition_order:
                self._condition_order.append(cond)

        return rtrn



class TargetedDesign(ExperimentalDesign):

    _samples: List[TargetedSample]
    _target_order: List[str]
    _treatment_order: List[str]

    def __init__(self, sample_to_condition: Dict[str, str],
                 sample_to_mark: Dict[str, str],
                 sample_to_treatment: Dict[str, str],
                 sample_to_rep_number: Dict[str, int] = None) -> None:
        super().__init__(sample_to_condition, sample_to_rep_number)

        # Convert Normal Samples to Targeted Samples
        targeted_samples = []
        for sample in self._samples:
            sample_name = sample.sample_name
            condition = sample.condition
            replicate = sample.replicate
            mark = sample_to_mark[sample_name]
            treatment = sample_to_treatment[sample_name]
            targeted_sample = TargetedSample(sample_name, condition, replicate,
                                             target=mark, treatment=treatment)
            targeted_samples.append(targeted_sample)

        self._target_order = []
        self._treatment_order = []

        self._samples = targeted_samples

    def get_samples(self, conditions: List[str] = None,
                    marks: List[str] = None,
                    treatments: List[str] = None,
                    reps: List[int] = None,
                    num_expected: int = None) -> List[str]:

        filters = args_to_filters(
            conditions=conditions,
            marks=marks,
            treatments=treatments,
            reps=reps
        )

        res = self.query('sample_name', filters)

        if num_expected and len(res) != num_expected:
            f_str = '\t\n' + '\t\n'.join(
                [s + ' - ' + self.get_condition(s)
                 for s in res])
            raise ValueError(f"Unexpected number of files found,"
                             f" {len(res)} != {num_expected} {f_str}")

        return res

    def set_target_order(self, target_order: List[str]) -> None:
        """
        Sets the order in which to sort the samples in this design.
        :param target_order: The order of the sample targets. Must 
            contain targets for all samples.
        :return:
        """
        check_one_to_one(target_order,
                         [sample.target for sample in self._samples],
                         'target')

        self._target_order = target_order

        for sample in self._samples:
            sample.set_target_order(target_order)

    def set_treatment_order(self, treatment_order: List[str]) -> None:
        """
        Sets the order in which to sort the samples in this design.
        :param treatment_order: The order of the sample treatments. Must 
            contain treatments for all samples.
        :return:
        """
        check_one_to_one(treatment_order,
                         [sample.treatment for sample in self._samples],
                         'treatment')
        self._treatment_order = treatment_order

        for sample in self._samples:
            sample.set_treatment_order(treatment_order)

    def get_mark(self, sample_name: str) -> str:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                return sample.target

    def get_target(self, sample_name: str) -> str:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                return sample.target

    def get_file_to_control(self, files: List[Path],
                            control_id: str = 'IgG') -> Dict[Path, Path]:
        file_to_control = {}
        for bam in files:
            sample_name = bam.name.split('_')[1]

            if self.get_mark(sample_name) == control_id:
                continue
            treatment = self.get_treatment(sample_name)
            control = self.find_in_files(files,
                                         treatments=[treatment],
                                         marks=[control_id],
                                         num_expected=1)[
                0]
            file_to_control[bam] = control
        return file_to_control

    def get_marks(self, samples: List[str] = None) -> List[str]:
        if not samples:
            samples = self.get_samples()
        return [self.get_mark(sample) for sample in samples]

    def get_treatment(self, sample_name: str) -> str:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                return sample.treatment

    def get_treatments(self, samples: List[str] = None) -> List[str]:
        if not samples:
            samples = self.get_samples()
        return [self.get_treatment(sample) for sample in samples]

    def __str__(self) -> str:
        def pretty_dict(dic: Dict[str, Any]) -> str:
            s = ''
            for sample in self._samples:
                s += '\t' + sample.sample_name + ': ' + str(
                    dic[sample.sample_name]) + '\n'
            return s

        return \
            f'=== Sample: Condition ===\n' \
            f'{pretty_dict({sample.sample_name: sample.condition for sample in self._samples})}\n' \
            f'=== Sample: Replicate ===\n' \
            f'{pretty_dict({sample.sample_name: sample.replicate for sample in self._samples})}\n' \
            f'=== Sample: Mark ===\n' \
            f'{pretty_dict({sample.sample_name: sample.target for sample in self._samples})}\n' \
            f'=== Sample: Treatment ===\n' \
            f'{pretty_dict({sample.sample_name: sample.treatment for sample in self._samples})}\n'

    def remove_sample(self, sample_name: str) -> None:
        for sample in self._samples:
            if sample.sample_name == sample_name:
                self._samples.remove(sample)
                return
        raise ValueError(f'{sample_name} not found in list of samples.')

    def __add__(self, other):
        """
        Returns an experimental design containing samples from both
        experimental designs. Resets condition orders if they are not
        identical.
        :param other: Another experimental design.
        :return:
        """

        if type(self) != type(other):
            raise ValueError('Addition not defined between these classes.')

        rtrn = deepcopy(self)

        rtrn._samples.extend(other._samples)

        for cond in other._target_order:
            if cond not in self._target_order:
                self._target_order.append(cond)

        for cond in other._treatment_order:
            if cond not in self._treatment_order:
                self._treatment_order.append(cond)

        return rtrn


def convert_to_targeted(design: ExperimentalDesign, mark_slice: slice,
                        treatment_slice: slice) -> TargetedDesign:
    """
    Converts the given experimental design into a targeted design.
    :param treatment_slice:
    :param mark_slice:
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

    return TargetedDesign(sample_to_cond,
                          sample_to_mark,
                          sample_to_treat,
                          sample_to_rep_num)


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

    for cond, rep in cond_rep_to_samples:
        samples = cond_rep_to_samples[(cond, rep)]

        # Split sample that have already been merged to maintain order.
        s_split = []
        for sample in samples:
            s_split.extend(sample.split('-'))
        samples = s_split

        samples.sort()
        sample = '-'.join(samples)

        sample_to_condition[sample] = cond
        sample_to_rep_number[sample] = rep

    return ExperimentalDesign(sample_to_condition,
                              sample_to_rep_number)


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
                              sample_to_rep_number)


def align_to_strs(to_align: Union[List[str], List[Path]],
                  strs: List[str]) \
        -> Union[List[Path], List[str]]:
    """
    Given that each sample_str contains a sample name, allowing for a 1 to 1
    mapping between each of to_align and strs, returns files ordered by
    strs.

    Does not consider contents of files past the first '.'.

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
                s_to_align = p_to_align.name.split('.')[0]
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
