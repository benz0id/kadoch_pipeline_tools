"""
Benjamin Tudor Price
10-02-2023

Command line access for <fix_sample_sheet>.

See <fix_sample_sheet> for detailed description of functionality.
"""
import argparse
from pathlib import Path
from components.sample_sheets import fix_sample_sheet
import logging
logging.root.disabled = True

parser = argparse.ArgumentParser("Correct Sample Sheet Formatting and "
                                 "get reverse compliemnts.")

parser.add_argument('sample_sheet', type=str,
                    help='Path to the sample sheet to be reformatted.')
parser.add_argument('-o', '--output-dir', type=str,
                    help='Path to the output output_dir. Parent of '
                         '<sample_sheet> by default.', default=None)
parser.add_argument('-f', '--rev-ifive', action='store_true',
                    help="Convert i5 indices to their reverse compliments.")
parser.add_argument('-s', '--rev-iseven', action='store_true',
                    help="Convert i7 indices to their reverse compliments.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Display additional information about progress.')

args = parser.parse_args()
if isinstance(args.output_dir, str):
    args.output_dir = Path(args.output_dir)

sheet = fix_sample_sheet(sample_sheet_path=Path(args.sample_sheet),
                         directory=args.output_dir,
                         rev_i5=args.rev_ifive,
                         rev_i7=args.rev_iseven,
                         verbose=args.verbose)

print(str(sheet))
exit(0)
