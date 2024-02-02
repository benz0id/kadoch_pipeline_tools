import itertools
import logging
import os
import subprocess
from copy import copy
from pathlib import Path
from typing import List, Callable, Tuple

import numpy as np
import qnorm as qnorm
from pandas import unique
from scipy.stats import stats

from utils.cache_manager import CacheManager
from utils.fetch_files import get_unique_filename, outpath_to_dirname
from sklearn.preprocessing import StandardScaler
from utils.job_formatter import ExecParams, PythonJob
from utils.job_manager import JobManager
from utils.path_manager import PathManager, cmdify
from constants.data_paths import HG19_IDXSTATS
from utils.utils import ExperimentalDesign
from operator import itemgetter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils.utils import ExperimentalDesign

HG19_GENOME = HG19_IDXSTATS

logger = logging.getLogger(__name__)

MAX_PCA_DIMS = 5


def generate_pca_plot(counts_matrix_path: Path,
                      design: ExperimentalDesign,
                      out_filepath: Path, dims: int = 3,
                      n_info_cols: int = 1, sample_ids: bool = False,
                      samples: List[str] = None,
                      colour_groups: List[str] = 'by_condition',
                      shape_groups: List[str] = 'by_rep',
                      title: str = 'PCA',
                      verbose: bool = False,
                      n_most_variable: int = None,
                      colour_label: str = 'labels',
                      shape_label: str = 'reps') -> sns.scatterplot:
    """
    Generates a PCA plot displaying a dimensionality reduced -
    representation of the given counts matrix.

    :param samples:
    :param shape_groups:
    :param colour_groups:
    :param dims:
    :param out_filepath:
    :param design: The design of the experiment.
    :param counts_matrix_path: A matrix with the number of counts for each
        peak.
    :param n_info_cols: The number of leading info columns to ignore.
    :param sample_ids: Whether to expect the column headers to be sample ids
        rather than sample names.
    :return:
    """
    # Extract raw data from the counts matrix.
    counts_dataframe = pd.read_csv(counts_matrix_path, sep='\t')

    matrix_samples = counts_dataframe.columns
    if n_info_cols > 0:
        to_rem = [matrix_samples[i] for i in range(n_info_cols)]
        counts_dataframe = counts_dataframe.drop(columns=to_rem)
        matrix_samples = counts_dataframe.columns
        counts_dataframe = counts_dataframe[matrix_samples].astype(float)

    if sample_ids:
        matrix_samples = [sample.strip().split('_')[1] for sample in matrix_samples]
        counts_dataframe.columns = matrix_samples

    if not samples:
        reps = [design.get_rep_num(label) for label in matrix_samples]
        conds = [design.get_condition(label) for label in matrix_samples]
        samples = matrix_samples
    else:
        reps = [design.get_rep_num(label) for label in samples]
        conds = [design.get_condition(label) for label in samples]
        counts_dataframe = counts_dataframe[samples]

    samples = list(samples)

    if isinstance(colour_groups, list):
        assert len(colour_groups) == len(samples)
    if isinstance(shape_groups, list):
        assert len(shape_groups) == len(samples)

    if n_most_variable:
        row_variability = np.std(counts_dataframe, axis=1)
        counts_dataframe = counts_dataframe.iloc[row_variability.argsort()[::-1][:n_most_variable]]

    counts_matrix = np.log2(counts_dataframe + 1)

    # norm_counts = scaler.fit_transform(counts_matrix)
    norm_counts = stats.zscore(counts_matrix, axis=0)

    pcs = qnorm.quantile_normalize(norm_counts.T)

    pca = PCA(min(MAX_PCA_DIMS, len(samples)))
    pcs = pca.fit_transform(pcs)

    pcdf = pd.DataFrame(data=pcs,
                        columns=['principal component ' + str(i)
                                 for i in range(1, pcs.shape[1] + 1)])
    pcdf['samples'] = samples

    if colour_groups == 'by_condition':
        pcdf['labels'] = conds
    if shape_groups == 'by_rep':
        pcdf['reps'] = reps

    if isinstance(colour_groups, list):
        pcdf['labels'] = colour_groups
    if isinstance(shape_groups, list):
        pcdf['reps'] = shape_groups

    pcdf.rename({'labels': colour_label,
                 'reps': shape_label}, inplace=True)

    if verbose:
        print(pcdf.iloc[:, -3:])

    logger.info((f'PCA Dataframe: for {counts_matrix_path}:' + '\n' +
                 str(pcdf)))

    props = pca.explained_variance_ratio_ * 100

    # Compare multiple principle components.
    for i, j in itertools.combinations(range(dims), 2):

        kwargs = {
            "data": pcdf,
            "x": 'principal component ' + str(i + 1),
            "y": 'principal component ' + str(j + 1)
        }

        if colour_groups is not None:
            kwargs['hue'] = 'labels'
            kwargs['palette'] = \
                sns.color_palette('bright', len(unique(pcdf['labels'])))
        if shape_groups is not None:
            kwargs["style"] = 'reps'

        sns.scatterplot(**kwargs)
        plt.xlabel(f'PC{i + 1}:{props[i]: .2f}%')
        plt.ylabel(f'PC{j + 1}:{props[j]: .2f}%')
        plt.title(title)

        filename = ''.join(['pc', str(i + 1), '_vs_', 'pc',
                            str(j + 1) + '.svg'])

        plt.savefig(out_filepath / filename, bbox_inches='tight')
        plt.close()






