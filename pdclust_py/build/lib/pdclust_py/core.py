import os
import itertools
import multiprocessing as mp

import numba
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time

@numba.jit(nopython=True)
def numba_pairwise_diss(np_arr1, np_arr2):
    print(np_arr1)

def read_bed_files(bed_files):

    csv_bed1 = pd.read_csv(bed_files[0], sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    csv_bed1['chromosome'] = csv_bed1['chr'].str.extract("([-+]?\d*\.\d+|[-+]?\d+)").astype(float)
    csv_bed1.drop(columns=['chr'])

    csv_bed2 = pd.read_csv(bed_files[1], sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    csv_bed2['chromosome'] = csv_bed2['chr'].str.extract("([-+]?\d*\.\d+|[-+]?\d+)").astype(float)
    csv_bed2.drop(columns=['chr'])

    if len(csv_bed1[csv_bed1['chromosome'].isna()]) > 1:
        print(bed_files[0])
    if len(csv_bed1[csv_bed1['chromosome'].isna()]) > 1:
        print(bed_files[1])

    # arr_bed1 = np.array(csv_bed1, dtype={'names': ['chr', 'pos', 'frac_methyl'], 'formats': [np.int16, np.int64, np.float16]})
    # arr_bed2 = np.array(csv_bed2, dtype={'names': ['chr', 'pos', 'frac_methyl'], 'formats': [np.int16, np.int64, np.float16]})
    # numba_pairwise_diss(arr_bed1, arr_bed2)
    # return arr_bed1, arr_bed2


def get_pairwise_dissimilarity(args):
    """
    Calculate pairwise dissimilarity between two cells
    :param args: tuple of two .bed files where each bed file contains the columns 'chromosome', 'start', 'end', 'methylation' (in percent, e.g. 0 or 100)
    :return: tuple of (row_label, col_label, pairwise_dissimilarity) where row_label is the name of the first cell, col_label is the name of the second cell, and pairwise_dissimilarity is the hamming distance between the two cells
    """

    row_label, col_label = args
    df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    row_label, col_label = row_label.split('/')[-1].replace('.bed', ''), col_label.split('/')[-1].replace('.bed', '')
    merged_df = pd.merge(df1, df2, on=['chr', 'pos'], suffixes=('_df1', '_df2'), how='inner')
    total_cpgs = merged_df.shape[0]
    pairwise_dissimilarity_total = np.sum(np.abs(merged_df['frac_methyl_df1'] - merged_df['frac_methyl_df2']))
    if total_cpgs != 0:
        pairwise_dissimilarity = pairwise_dissimilarity_total / total_cpgs
    else:
        pairwise_dissimilarity = 0.0
    return row_label, col_label, pairwise_dissimilarity


def make_unique_combinations(path_to_bed_files):
    """
    Make unique combinations of all cells in a directory
    :param path_to_bed_files: containing the BED files
    :return: list of tuples of unique combinations of all cells
    """
    if not path_to_bed_files.endswith('/'):
        path_to_bed_files = path_to_bed_files + '/'
    file_list = [path_to_bed_files + f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
    assert len(file_list) != 0, f'No .bed files found in "{path_to_bed_files}"'
    unique_combinations = []
    for file_pair in itertools.combinations(file_list, 2):
        if set(file_pair) not in [set(pair) for pair in unique_combinations]:
            unique_combinations.append(file_pair)
    print('Generation of unique combinations of cells done')
    return unique_combinations


# def pdclust(path_to_bed_files, plot=False, heatmap_output_path='heatmap.png'):
#     """
#     Calculate pairwise dissimilarity between all cells (in BED file format) in a directory
#     :param path_to_bed_files:  containing the BED files
#     :return: pandas dataframe containing pairwise dissimilarity between all cells
#     """
#     unique_combinations = make_unique_combinations(path_to_bed_files)
#     columns = [x.split('/')[-1].replace('.bed', '') for (x,y) in unique_combinations]
#     index = [x.split('/')[-1].replace('.bed', '') for (x,y) in unique_combinations]
#     result_df = pd.DataFrame(index=index, columns=columns)
#     with mp.Pool(mp.cpu_count()) as pool:
#         results = pool.map(get_pairwise_dissimilarity, unique_combinations)
#     for (row, col, pair_dist) in results:
#         result_df.at[row, col] = pair_dist
#         result_df.at[col, row] = pair_dist
#         result_df.at[col, col] = 0.0
#     result_df.fillna(0.0, inplace=True)
#     if plot: plot_heatmap(result_df)
#     return result_df


def trivial_pdclust(path_to_bed_files, plot=True):
    if not path_to_bed_files.endswith('/'):
        path_to_bed_files = path_to_bed_files + '/'

    file_list = [path_to_bed_files + f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
    pool = mp.Pool(processes=mp.cpu_count())
    # pool = mp.Pool(processes=4)

    results = pool.map(get_pairwise_dissimilarity,
                       [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [x.split('/')[-1].replace('.bed', '') for x in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for i, elem in enumerate(results):
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity
    result_df.fillna(0.0, inplace=True)
    if plot:
        plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap.png')
    pool.close()
    return result_df


def plot_heatmap(dist_matrix, output_path='./heatmap.png'):
    """
    Plot heatmap of pairwise dissimilarity between all cells with annotation bars depending on cell types
    :param output_path: path to save the heatmap
    :param dist_matrix: return value from pdclust(). Alternatively: Distance matrix between all cells (pandas dataframe)
    :return: None, saves heatmap at the destination of output_path
    """
    # colors = ['red' if 'ESLAM' in x else 'blue' for x in dist_matrix.columns]
    cluster_map = sns.clustermap(data=dist_matrix, method='ward', metric='euclidean', figsize=(20, 20), cmap='YlOrRd',
                                 xticklabels=False, yticklabels=False)
    # vmin=12.5, vmax=17)

    cluster_map.ax_row_dendrogram.set_visible(False)
    cluster_map.ax_col_dendrogram.set_visible(True)

    plt.savefig(output_path)


if __name__ == '__main__':
    trivial_pdclust('/home/ptheo/Documents/university/bachelor/data/human_data/')
