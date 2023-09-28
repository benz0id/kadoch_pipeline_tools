### kmeans cluster helper fx that
### returns sorted clustered df with
### cluster labels
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from seaborn.matrix import ClusterGrid
from sklearn.cluster import KMeans
import seaborn as sns
import numpy as np


def load_counts_matrix(path: Path, sep: str = '\t') -> pd.DataFrame:
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
    plt.savefig(outfile)
    return g


def write_out_bed(df: pd.DataFrame, out_path: Path) -> None:
    """
    Writes out a dataframe as a bedfile.
    :param df: A dataframe, with the row names as bedfiles.
    :param out_path: Path to the output bedfile.
    :return: None
    """
    # TODO something like this
    with open('atac.clusters.bed', 'w') as outfile:
        for i in df_k_reordered.index:
            coordinate = i.strip().replace(':', '\t').replace('-', '\t')
            # we want to change chr1:500-1000 to chr1\t500\t1000
            outfile.write(coordinate + '\n')



