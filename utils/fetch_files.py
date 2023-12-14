import logging
import os
import sys
import re
import uuid
from copy import copy
from datetime import datetime
from pathlib import Path
from typing import List, Union

from utils.path_manager import cmdify, PathManager

logger = logging.getLogger(__name__)

Regex = str
Regexes = Union[List[Regex], Regex]
def get_matching_strs(strs: List[str],
                      matching: Regexes,
                      not_matching: Regexes = None,
                      containing: bool = False,
                      under_delim: bool = False,
                      verbose: bool = False,
                      filetype: str = None,
                      n: int = None,
                      one_to_one: bool = False,
                      inds: bool = False) -> Union[List[str], List[int]]:
    """
    Gets all strings in <strs> that match any of <matching> and none of
    <not_matching>.
    :param strs: A list of strings.
    :param matching: String(s) representing regular expressions.
    :param not_matching: String(s) representing regular expressions.
    :param containing: Check whether a match to regexes exists in a substring
     of each of <strs>, rather than an exact match.
    :param under_delim: All matching strings must be wrapped contained between '_',
       or at the start or end of the string. For use with filenames with elements
       seperated by '_'.
    :param verbose: Whether to print information regarding the process to the
     console.
    :param inds: Whether to return the indecies of the matches rather than the
        actual matches themselves.
    :param n: Fails unless exactly <n> matches are found.
    :return: All strings in <strs> matching <matching> and not matching
    <not_matching>.
    """
    # Do some typechecking to allow for iterable or non-iterable inputs.
    if not_matching is None:
        not_matching = []
    elif not isinstance(not_matching, list):
        not_matching = [not_matching]
    if not isinstance(matching, list):
        matching = [matching]
    
    # Only needed if one-to-one mapping is required.
    inds_matches = [False for _ in matching]

    if under_delim:
        matching = ['(^|_)' + regex + '(_|$)' for regex in matching]
        not_matching = ['.*' + regex + '.*' for regex in not_matching]

    # Change regexes to search for matching substring.
    if containing or under_delim:
        matching = ['.*' + regex + '.*' for regex in matching]
        not_matching = ['.*' + regex + '.*' for regex in not_matching]

    if filetype:
        for i in range(len(matching)):
            matching[i] = matching[i] + '.*' + filetype + '$'

    valid = []
    # Find all matching strings.
    for s in strs:
        matches = [re.match(reg, s) for reg in matching]
        matches_wanted_re = any(matches)

        matches = [re.match(reg, s) for reg in not_matching]
        matches_not_wanted_re = any(matches)

        if not matches_wanted_re or matches_not_wanted_re:
            continue

        valid.append(s)

        # Ensure that the added string did not match multiple regexes.
        print(matches)
        print(sum(matches))
        if one_to_one and sum(matches) != 1:
            RuntimeError(
                f"One-to-one mapping not found. Multiple matches found "
                f"for {s}.")

        # Ensure that the added string does not match any of the previously
        # matched regexes.
        for i, match in enumerate(matches):
            if one_to_one and match and inds_matches[i]:
                raise RuntimeError(f"One-to-one mapping not found. Multiple "
                                   f"matches found for {matching[i]}.")

        inds_matches = [inds_matches[i] or matches[i]
                        for i in range(len(inds_matches))]

    s = ' '.join([
        '\nSearching for matches in \n\t', '\n\t'.join(strs),
        '\nMust match one of: ', ', '.join(matching),
        '\nCannot match any of: ', ', '.join(not_matching),
        '\nMatches found: \n\t', '\n\t'.join(valid)
    ])
    logger.debug(s)

    if verbose:
        print(s)

    if inds:
        vald_inds = []
        strs_copy = copy(strs)
        for v in valid:
            ind = strs_copy.index(v)
            strs_copy[ind] = ''
            vald_inds.append(ind)
        valid = vald_inds

    if n is not None and len(valid) != n:
        logger.debug(s)
        print(s, file=sys.stderr)
        raise RuntimeError(f'Unexpected number of files found. \n\tExpected: {n} \n\tFound: {len(valid)}')

    return valid


def get_matching_files(directory: Path,
                       matching: Regexes,
                       not_matching: Regexes = None,
                       filetype: str = None,
                       containing: bool = False,
                       under_delim: bool = False,
                       paths: bool = False,
                       n: int = None,
                       one_to_one: bool = False,
                       verbose: bool = False) -> Union[List[str], List[Path]]:
    """
    Gets all filenames in the given <directory> that match any of <matching>
    and none of <not_matching>.
    :param directory: The directory to be searched.
    :param matching: String(s) representing regular expressions.
    :param paths: Whether to convert all
    :param not_matching: String(s) representing regular expressions.
    :param containing: Check whether a match to regexes exists in a substring
     of each of <strs>, rather than an exact match.
    :param under_delim: All matching strings must be wrapped contained between '_',
        or at the start or end of the string. For use with filenames with elements
        seperated by '_'.
    :param verbose: Whether to print information regarding the process to the
     console.
    :param n: Fails unless exactly <n> matches are found.
    :return: All strings in <strs> matching <matching> and not matching
    <not_matching>.
    """
    if not directory.exists():
        raise ValueError(str(directory) + " does not exist.")
    if not directory.is_dir():
        raise ValueError(str(directory) + " is not a directory.")

    matching_files = get_matching_strs(os.listdir(directory),
                                       matching,
                                       not_matching,
                                       containing=containing,
                                       under_delim=under_delim,
                                       verbose=verbose,
                                       filetype=filetype,
                                       n=n,
                                       one_to_one=one_to_one)

    if paths:
        matching_files = [Path(directory / match)
                          for match in matching_files]

    return matching_files


def copy_to(directory: Path, file: Union[Path, List[Path]],
            avoid_recopy: bool = False) -> Path:
    """
    Copies the given file into the given directory. Overwrites any files
    with the same name in the directory.
    :param avoid_recopy: Do no copy files that exist at target dest iff true.
    :param file: A file.
    :param directory: A directory.
    :return: Filepath to the new file.
    """
    if isinstance(file, list):
        for ind_file in file:
            copy_to(directory,
                    ind_file,
                    avoid_recopy)

    new_path = directory / file.name

    if avoid_recopy and new_path.exists():
        return new_path

    os.system(' '.join(['cp', str(file), str(new_path)]))
    return new_path


def copy_to_cmds(directory: Path, files: List[Path],
                 avoid_recopy: bool = False) \
        -> List[str]:
    """
    Copies the given file into the given directory. Overwrites any files
    with the same name in the directory.
    :param avoid_recopy: Do no copy files that exist at target dest iff true.
    :param file: A file.
    :param directory: A directory.
    :return: Filepath to the new file.
    """
    cmds = []
    for file in files:
        if (directory / file.name).exists() and avoid_recopy:
            pass
        else:
            cmds.append(cmdify('cp', file, directory))
    return cmds


def get_unique_filename() -> str:
    """
    Gets a unique filename using the current time.
    :return: A unique filename.
    """
    return str(uuid.uuid4())


def outpath_to_dirname(path: Path) -> str:
    """
    Create a name for an output directory based on some filepath.
    :param path: Path to a file.
    :return: A unique directory name based on the input <path>
    """
    out = str(path).replace('/', '.')

    if len(out) > 255:
        raise ValueError("Directory names cannot exceed 255 characters.")
    return out


def extract_figs(org_dir: Path, new_dir: Path, path_manager: PathManager,
                 filetypes: List[str] = None, cmds: List[str] = None,
                 depth: int = 0, verbose: bool = False) -> List[str]:
    """
    Copy all figures of <filetypes> from <org_dir> to a mirror directory in <new_dir>.
    :param org_dir: The original directory.
    :param new_dir: The new directory to be created. Will have the exact same
        file structure as org_dir, but containing only files with one of
        <filetypes>.
    :param path_manager: Path manager to aid with directory creation.
    :param filetypes: The filetypes to copy.
    :param cmds: Used for passing commands recursively. Can be ignored.
    :return:
    """
    if not filetypes:
        filetypes = ['pdf', 'svg', 'tiff', 'png', 'jpg', 'jpeg']
    if cmds is None:
        cmds = []

    if verbose:
                print('\t' * depth, org_dir.name)

    for f in os.listdir(org_dir):
        sub_org = org_dir / f
        ext = f.split('.')[-1]

        if sub_org.is_dir():
            sub_new = new_dir / f
            extract_figs(sub_org, sub_new, path_manager, filetypes, cmds, depth + 1, verbose)
        elif ext in filetypes:
            if verbose:
                print('\t' * depth, sub_org.name)
            new_dir = path_manager.make(new_dir)
            cmds.extend(copy_to_cmds(new_dir, [sub_org]))

    return cmds


