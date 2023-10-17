import logging
import os
import re
import uuid
from copy import copy
from datetime import datetime
from pathlib import Path
from typing import List, Union

from utils.path_manager import cmdify

logger = logging.getLogger(__name__)

Regex = str
Regexes = Union[List[Regex], Regex]
def get_matching_strs(strs: List[str],
                      matching: Regexes,
                      not_matching: Regexes = None,
                      containing: bool = False,
                      verbose: bool = False,
                      inds: bool = False) -> Union[List[str], List[int]]:
    """
    Gets all strings in <strs> that match any of <matching> and none of
    <not_matching>.
    :param strs: A list of strings.
    :param matching: String(s) representing regular expressions.
    :param not_matching: String(s) representing regular expressions.
    :param containing: Check whether a match to regexes exists in a substring
     of each of <strs>, rather than an exact match.
    :param verbose: Whether to print information regarding the process to the
     console.
    :param inds: Whether to return the indecies of the matches rather than the
        actual matches themselves.
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

    # Change regexes to search for matching substring.
    if containing:
        matching = ['.*' + regex + '.*' for regex in matching]
        not_matching = ['.*' + regex + '.*' for regex in not_matching]

    valid = []
    # Find all matching strings.
    for s in strs:
        matches = [re.match(reg, s) for reg in matching]
        matches_wanted_re = any(matches)

        matches = [re.match(reg, s) for reg in not_matching]
        matches_not_wanted_re = any(matches)

        if matches_wanted_re and not matches_not_wanted_re:
            valid.append(s)

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

    return valid


def get_matching_files(directory: Path,
                       matching: Regexes,
                       not_matching: Regexes = None,
                       containing: bool = False,
                       paths: bool = False,
                       verbose: bool = False) -> Union[List[str], List[Path]]:
    """
    Gets all filenames in the given <directory> that match any of <matching>
    and none of <not_matching>.
    :param directory: The directory to be searched.
    :param matching: String(s) representing regular expressions.
    :param not_matching: String(s) representing regular expressions.
    :param paths: Whether to convert all
    :param containing: Check whether a match to regexes exists in a substring
     of each of <strs>, rather than an exact match.
    :param verbose: Whether to print information regarding the process to the
     console.
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
                                       containing,
                                       verbose)

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


