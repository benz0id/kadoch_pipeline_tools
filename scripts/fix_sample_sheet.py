"""
Benjamin Tudor Price
10-02-2023

See fix_sample_sheet for detailed description of functionality.
"""
import argparse
from pathlib import Path
from components.sample_sheets import fix_sample_sheet
import logging
logging.root.disabled = True


parser = argparse.ArgumentParser("Correct Sample Sheet Formatting and get "
                                 "Index2 Reverse Complements")

parser.add_argument('sample_sheet', type=str,
                    help='Path to the sample sheet to be reformatted.')
parser.add_argument('-o', '--output-dir', type=str,
                    help='Path to the output output_dir. Parent of '
                         '<sample_sheet> by default.', default=None)
parser.add_argument('-n', '--no-rev', action='store_true',
                    help="Do not apply reverse compliment to i5 indexes.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Display additional information about progress.')

args = parser.parse_args()
if isinstance(args.output_dir, str):
    args.output_dir = Path(args.output_dir)

sheet = fix_sample_sheet(sample_sheet_path=Path(args.sample_sheet),
                         directory=args.output_dir,
                         rev_comp=not args.no_rev,
                         verbose=args.verbose)

print(str(sheet))
