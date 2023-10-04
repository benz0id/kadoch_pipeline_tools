import logging
from copy import copy
from pathlib import Path
from typing import List, Dict, Callable

from utils.fetch_files import get_unique_filename
from utils.job_formatter import ExecParams
from utils.job_manager import JobManager
from utils.path_manager import cmdify


logger = logging.getLogger(__name__)


class MultiIntersector:
    """

    === Description ===
    Get the intersections from an arbitrary number of bedfiles.


    === Private attributes ===

    bedfiles
        The bedfiles to be intersected.
    bed_identifiers
        Names to write out when creating intersection files.
    binary_counter
        A binary counter with length equal to that of <bedfiles>.
    temp_directory
        Directory to store temporary files.
    """

    _bedfiles: List[Path]
    _bed_identifiers: List[str]
    _binary_counter: List[bool]
    _temp_directory: Path
    _light_job: ExecParams
    _jobs: JobManager
    _run: Callable

    def __init__(self, temp_directory: Path,
                 job_manager: JobManager,
                 light_job: ExecParams) -> None:
        self._bedfiles = []
        self._bed_identifiers: List[str]
        self._binary_counter: List[bool]

        self._jobs = job_manager
        self._light_job = light_job
        self._temp_directory = temp_directory

        self._run = lambda x: self._jobs.execute(x, self._light_job)

    def iter_counter(self) -> bool:
        """
        Iterates the binary counter.
        :return: Whether the counter is complete.
        """

        i = 0
        while self._binary_counter[i]:
            self._binary_counter[i] = False
            i += 1

        if i == len(self._binary_counter):
            return False
        else:
            self._binary_counter[i] = True
            return True

    def get_sorted_intersection(self, beds: List[Path]) -> Path:
        """
        Get the peaks in beds that overlap at least one peak in each other
        bedfile. Unmerged, Sorted.

        :param beds: A list of bedfiles.
        :return: Path to sorted bedfile.
        """
        inters = []
        set_name = '_'.join([b.name.split('_')[1] for b in beds])

        # Do not need to intersect beds if we only have one.
        if len(beds) == 1:
            return beds[0]

        for bed in beds:
            inter_name = bed.name[:-4] + '-inter-' + set_name + '.bed'

            inter = self._temp_directory / inter_name
            inters.append(inters)

            others = copy(beds)
            others.remove(bed)
            print(beds)

            cmd = cmdify(
                'bedtools intersect -a', bed,
                '-b', *others,
                '-wa > ', inter
            )
            self._run(cmd)

        final_name = set_name + '_sorted_intersection.bed'
        final_file = self._temp_directory / final_name
        cmd = cmdify(
            'cat', *inters,
            '| bedtools sort >', final_file
        )
        self._run(cmd)

        return final_file

    def get_all_merged_and_sorted(self, beds: List[Path]) -> Path:
        """
        Get all peaks in the given bedfiles after merging and sorting.
        :param beds: A list of bedfiles.
        :return: Path to sorted and merged bedfile.
        """
        set_name = get_unique_filename()
        out_path = self._temp_directory / (set_name + '_merged_and_sorted.bed')

        cmd = cmdify(
            'cat', *beds,
            '| bedtools sort',
            '| bedtools merge',
            '>', out_path
        )
        self._run(cmd)

        return out_path

    def get_name(self) -> str:
        """
        Returns the name of the bedfile that would be created by generating
        a venn section from the current state of <self._binary_counter>
        :return: Name of a bedfile.
        """
        return '_'.join(self.get_identifiers())

    def get_identifiers(self) -> List[str]:
        """
        Get the identifiers corresponding to the currently selected bedfiles,
        as denoted by self._binary_counter.
        :return: Name of a bedfile.
        """
        included = []
        for i, include in enumerate(self._binary_counter):
            if include:
                included.append(self._bed_identifiers[i])
        return included

    def get_included_bedfiles(self) -> List[Path]:
        """
        Gets all bedfiles currently included in the given intersection, as
        denoted by self._binary_counter.
        :return: A list of included bedfiles.
        """
        included = []
        for i, include in enumerate(self._binary_counter):
            if include:
                included.append(self._bedfiles[i])
        return included

    def get_excluded_bedfiles(self) -> List[Path]:
        """
        Gets all bedfiles NOT currently included in the given intersection, as
        denoted by self._binary_counter.
        :return: A list of not included bedfiles.
        """
        not_included = []
        for i, include in enumerate(self._binary_counter):
            if not include:
                not_included.append(self._bedfiles[i])
        return not_included

    def update_counts_dict(self, bedfile: Path, counts_dict: Dict[str, int]
                           ) -> None:
        """
        Add the number of counts in the given file to the counts dict using the
        currently stored intersection as a key.
        :param bedfile: Path to a valid bedfile.
        :param counts_dict: A dictionary mapping the name of an intersection to
            the number of counts in that interseciton.
        """
        with open(bedfile, 'r') as inf:
            num_lines = len(inf.readlines())
        s = ' '.join(self.get_identifiers()) + '.bed'
        counts_dict[s] = num_lines

    def get_all_venn_intersections(self, bedfiles: List[Path],
                                   identifiers: List[str],
                                   output_dir: Path,
                                   merge: bool = True) -> Dict[str, int]:
        """
        Gets the possible exclusive set combinations of the given
        bedfiles, essentially giving all section of a venn diagram as peaks
        files.

        === Detailed Implementation Description ===

        Each section of the venn diagram is formed by the exclusive
        intersection (EI) of some number of bed files.

        Let us consider the intersection between bedfiles {A, B, C} from a
        venn diagram of {A, B, C, D, E}.
        In order for a peak to lie within this intersection it must

        1. Be contained within one of {A, B, C} (trivially).
        2. Overlap with at least one peak in each of {A, B, C}.
        3. Not overlap with any peak in any of {D, E}.

        All peaks matching these conditions are found, and then merged
        to produce the final bedfile. This is repeated for each combination
        of {A, B, C, D, E}

        :param bedfiles: A list of sorted bedfiles.
        :param identifiers: A list of strings, where each string is an
            identifier for the bedfile at the corresponding position
            in bedfiles.
        :param output_dir: The directory in which to place the merged
            bedfiles.
        :param merge: Whether to merge the final output bedfiles.
        :return: A dictionary mapping identifiers (seperated by spaces)
            to the number of counts in that exclusive intersection.
        """
        # Initialise attributes.
        self._binary_counter = [False] * len(bedfiles)
        self._bed_identifiers = copy(identifiers)
        self._bedfiles = copy(bedfiles)

        # Specify addition merge step if neccesary.
        if merge:
            merge_str = ' | bedtools merge'
        else:
            merge_str = ''

        counts_dict = {}

        # Generate exclusive intersections.
        while self.iter_counter():
            out_path = output_dir / (self.get_name() + '.bed')

            included = self.get_included_bedfiles()
            excluded = self.get_excluded_bedfiles()

            common_included = self.get_sorted_intersection(included)
            merged_excluded = self.get_all_merged_and_sorted(excluded)

            cmd = cmdify(
                'bedtools intersect',
                '-a', common_included,
                '-b', merged_excluded,
                '-wa',
                '| bedtools sort' + merge_str,
                '>', out_path
            )
            self._run(cmd)

            self.update_counts_dict(out_path, counts_dict)

        return counts_dict










