from pathlib import Path
from typing import List

from components.peak_counter import PeakCounter
from utils.utils import ExperimentalDesign


class BedSorter:
    """
    Responsible for the sorting of bedfiles.
    """

    _design: ExperimentalDesign
    _counts_dir: Path
    _peak_counter: PeakCounter

    def __init__(self, design: ExperimentalDesign, counts_dir: Path,
                 peak_counter: PeakCounter) -> None:
        """
        Create a bed sorter using the given attributes.
        :param design: The experimental design used to sort samples and
            find files matching sample names.
        :param counts_dir: A directory containing reads or fragments in the
            bam or bed formats (respectively). Must not contain both bam and
            bed files.
        """
        self._design = design
        self._counts_dir = counts_dir
        self._peak_counter = peak_counter

    def sort_bed(self, bed: Path, samples: List[str], out_path: Path) -> Path:
        """
        Sort a bedfile by the mean number of counts at each peak across all
        of <counts>.
        :param bed: A bedfile.
        :param samples: The samples of the reads files to use.
        :param out_path: The path to the output sorted bed file.
        :return: out_path
        """
        matrix_out = out_path.parent / ('.'.join(out_path.name.split('.')[:-1])
                                        + '.tsv')

        self._peak_counter.get_matrix(bed, self._counts_dir, self._design,
                                      matrix_out, samples)





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

