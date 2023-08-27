from pathlib import Path
import os
from typing import List

pd = Path(__file__).parent
bams_path = pd / 'bam'

CUTOFF = 0.01


def run(cmd: str) -> None:
    print(cmd)
    os.system(cmd)


def replace_mkdir(path: Path) -> None:
    if path.exists():
        run('rm -r ' + str(path))
    run('mkdir ' + str(path))


def get_all_matching(strs, matches,
                     not_matches='THIS WILL NEVER BE IN ANYTHING') -> List[
    str]:
    """Returns all of strs that match <matches> and not <not_matches>"""
    # print("\n\nGetting all matching - \n" + str(strs) + '\n' + matches + '\n' + not_matches)
    rtrn = []
    for s in strs:
        if matches in s and not_matches not in s:
            rtrn.append(s)
    # print("Got: " + str(rtrn) + '\n\n')
    return rtrn


def copy_to(filepaths: List[Path], dir):
    """Copes of of the given filepaths to the given directory."""
    for filepath in filepaths:
        file_name = str(filepath).split('/')[-1]
        run(' '.join([
            'cp', str(filepath), str(dir / file_name)
        ]))


treatments = ["G401", "NCIH1048_DMSO", "NCIH1048_FHD286"]
targets = ["CTCF", "H3K27ac", "POU2F3", "RAD21", "BRD9", "RPB1", "BRG1"]


def generate_peaks():
    inputs = get_all_matching(os.listdir(bams_path), 'Input')

    for treat in treatments:
        # Create useful directories.
        treat_bams_path = pd / (treat + '_bam')
        treat_peaks_path = pd / (treat + '_peaks')
        # replace_mkdir(treat_bams_path)
        # replace_mkdir(treat_peaks_path)

        # Copy over everything except input bams.
        # treat_bams = get_all_matching(os.listdir(bams_path), treat, 'Input')
        # copy_to([bams_path / treat_bam for treat_bam in treat_bams], treat_bams_path)

        # Get input bamfile.
        print("Inputs: ", inputs, '\n\nTreatments: ', treat_bams, '\n\n')
        input_bam = list(set(inputs).intersection(
            set(get_all_matching(os.listdir(bams_path), treat))))[0]
        print("Input: ", input_bam)

        # Run peak calling.
        run(' '.join([
            'python $soft_kc/Run_MACS2peak_calling.py',
            #            "-i", str(bams_path /input_bam),
            "-seq", 'atac',
            "-bams", str(treat_bams_path),
            "-o", str(treat_peaks_path),
            "-p", 'narrow',
            "-c", str(CUTOFF),
            "-g", "hs"
        ]))


def get_peak_counts(peaks_dir) -> None:
    for peak_file in os.listdir(peaks_dir):
        if 'narrowPeak' not in peak_file:
            continue

        with open(peaks_dir / peak_file, 'r') as pf:
            file_comps = peak_file.split('_')
            cond, pd = file_comps[4], file_comps[5]
            s = cond + ' ' + pd
            print(f'{s}', ' - ', len(pf.readlines()))


def generate_venns() -> None:
    DMSO_peaks_path = pd / "NCIH1048_DMSO_ht_peaks"
    FHD286_peaks_path = pd / "NCIH1048_FHD286_ht_peaks"

    for target in targets:
        pass


DMSO_ht_peaks_path = pd / "NCIH1048_DMSO_ht_peaks"
FHD286_ht_peaks_path = pd / "NCIH1048_FHD286_ht_peaks"

DMSO_lt_peaks_path = pd / "NCIH1048_DMSO_lt_peaks"
FHD286_lt_peaks_path = pd / "NCIH1048_FHD286_lt_peaks"


def main():
    print('\n\n=== DMSO Peaks 0.05 ===')
    get_peak_counts(DMSO_ht_peaks_path)
    print('\n\n=== FHD 286 Peaks 0.05 ===')
    get_peak_counts(FHD286_ht_peaks_path)

    print('\n\n=== DMSO Peaks 0.01 ===')
    get_peak_counts(DMSO_lt_peaks_path)
    print('\n\n=== FHD 286 Peaks 0.01 ===')
    get_peak_counts(FHD286_lt_peaks_path)


if __name__ == '__main__':
    main()

# Run_MACS2peak_calling.py [-h] [-i INPUT] [-seq {chip,atac}] -bams LISTOFBAMS -o DIRECTORY -p {broad,narrow} [-c CUTOFF] [-g {hs,mm,vero}]



