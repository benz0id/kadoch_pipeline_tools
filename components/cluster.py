### kmeans cluster helper fx that
### returns sorted clustered df with
### cluster labels
from copy import copy
from pathlib import Path
from typing import Union, List

import matplotlib.pyplot as plt
import pandas as pd
from seaborn.matrix import ClusterGrid
from sklearn.cluster import KMeans
from scipy import stats
import seaborn as sns
import numpy as np

from utils.utils import ExperimentalDesign


def load_counts_matrix(path: Path, sep: str = '\t',
                       design: ExperimentalDesign = None) -> pd.DataFrame:
    """
    For some counts-like matrix, with a leading column consisting of
    identifiers such as gene names or peaks, load that matrix as a pandas
    dataframe with those identifiers as the row names.
    :param path: The path to the counts-like file.
    :param sep: The seperator in the file.
    :return: A formatted dataframe.
    """
    df = pd.read_csv(path, sep=sep)
    df = df.set_index(df.columns[0])

    if design:
        # Run a bit of filtering and checking on expression data.
        invalid_columns = design.get_invalid_sample_inds(
            df.columns)
        to_remove = [df.columns[i] for i in invalid_columns]
        df = df.drop(columns=to_remove)

        reordered_columns = design.align_to_samples(df.columns)
        df = df[[reordered_columns]]
        df.columns = design.get_sample_descs()

    return df


def add_clustering(df: pd.DataFrame, n_clusters, random_state=42) -> \
        pd.DataFrame:
    """
    Add a column "cluster" to the dataframe, as determined by the
    sklearn.cluster.Kmeans algorithm
    :param df: A dataframe of numeric values.
    :param n_clusters: The number of clusters to assign.
    :param random_state: The random state to initialise clusters and shuffle
    data.
    :return:
    """
    # Shuffle dataframe.
    df = df.sample(frac=1, random_state=random_state)

    # Run K means.
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(df)

    # Add cluster column.
    return df.assign(clusters=kmeans.labels_).sort_values('clusters')


def vis_clustered_data(df: pd.DataFrame, outfile: Path,
                       **kwargs) -> ClusterGrid:
    """
    Create and save a heatmap visualising the given dataframe.

    :param df: A dataframe with a number of numeric columns, and a single
    "clusters" column (see add_cluster_col).

    :param outfile: The filepath to the output heatmap.

    :param kwargs: Some number of arguments to be passed to the clustermap
        function, overwrites funtion stored defaults.

    :return: The clustergrid heatmap object.
    """

    defaults = {
        "cmap": 'bwr',
        "center": 0,
        "vmin": -2,
        "vmax": 2,
        "row_cluster": False,
        "col_cluster": False,
        "row_colors": ['C%d' % i for i in df.clusters]
    }

    defaults.update(kwargs)
    g = sns.clustermap(df.drop(columns='clusters'), **defaults)
    plt.savefig(outfile, dpi=500)
    return g


def quick_clustering_analysis(expression_data: Union[Path, pd.DataFrame],
                              n_clusters: int,
                              out_path: Path,
                              design: ExperimentalDesign,
                              transform: str = 'zscore',
                              save_out: bool = False) -> pd.DataFrame:
    """
    Run clustering analysis on the given expression data.

    :param expression_data: Path to counts-like expression_data
        or pre-parsed counts-like data.

    :param n_clusters: The number of clusters to form.

    :param out_path: The path to the output heatmap image.

    :param design: The experimental design.

    :param transform: How to transform each row/value before running
        clustering. One of ["log", "zscore", "none"].

    :param save_out: Whether to save the clustered data as a matrix
        (clustered_rpkms.tsv) to <out_path>.

    :return: Dataframe with clustered reads.
    """

    if isinstance(expression_data, Path):
        expression_data = load_counts_matrix(expression_data, design=design)
    else:
        expression_data = copy(expression_data)

    if transform == 'zscore':
        expression_data = stats.zscore(expression_data, axis=1)
    elif transform == 'log':
        expression_data = np.log2(expression_data + 1)
    elif transform == 'none':
        pass
    else:
        raise ValueError("Unknown transform requested")

    clustered_data = add_clustering(expression_data,
                                    n_clusters=n_clusters)

    if save_out:
        clustered_data.to_csv(out_path / 'clustered_rpkms.tsv', sep='\t')

    vis_clustered_data(clustered_data, out_path)

    return clustered_data


def get_clusters(clustered_counts: pd.DataFrame,
                 cluster_nums: Union[int, List[int]],
                 out_dir: Path) -> Union[List[Path], Path]:
    """
    Extract all of <cluster_nums> from <clustered_counts>, format them as
    bedfiles and place them in <out_dir>.

    :param clustered_counts: A dataframe with an index formatted as a list of
        bed sites formatted <chr>:<start>-<stop> and a column named "clusters",
        containing the integer assignment to which each row belongs.

    :param cluster_nums: The numbers of the clusters to be extracted.
    :param out_dir: The directory into which the extracted bedfiles should be
        placed.

    
    """
    if isinstance(cluster_nums, int):
        cluster_nums = [cluster_nums]
        unpack = True
    else:
        unpack = False


    out_peak_files = []
    for cluster_num in cluster_nums:
        out_file = out_dir / f'peak_number_{cluster_num}.bed'
        out_peak_files.append(out_file)
        out = open(out_file, 'w')

        for index in clustered_counts.index:
            cluster = clustered_counts[index, 'clusters']
            if cluster == cluster_num:
                coordinate = index.\
                    strip().\
                    replace(':', '\t').\
                    replace('-', '\t')
                out.write(coordinate + '\n')
        out.close()

    if unpack:
        return out_peak_files[0]
    return out_peak_files





