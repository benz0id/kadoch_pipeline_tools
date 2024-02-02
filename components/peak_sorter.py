from pathlib import Path
from typing import List


class BedSorter:
    """
    Responsible for the sorting of bedfiles.
    """

    def sort_bed(self, bed: Path, counts: List[Path], out_path: Path) -> Path:
        """
        Sort a bedfile by the mean number of counts at each peak across all
        of <counts>.
        :param bed: A bedfile.
        :param counts: A list of bed/bam files.
        :param out_path: The path to the output sorted bed file.
        :return: out_path
        """
        pass

    def comparison_sort(self, bed: Path, count1: Path, counts2: Path,
                        out_path: Path, metric: str = 'fold_change') -> Path:
        """
        Sorts <bed> according to some feature of
        :param bed:
        :param count1:
        :param counts2:
        :param out_path:
        :param metric:
        :return:
        """

