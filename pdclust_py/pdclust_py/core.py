import os
import itertools
import multiprocessing as mp

import numba
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# @numba.jit(forceobj=True)
# def numba_pair_dist(chr_arr1, chr_arr2, pos_arr1, pos_arr2, meth_arr1, meth_arr2):
#     result = count = idx1 = idx2 = 0
#
#     while idx1 < len(chr_arr1):
#         while idx2 < len(chr_arr2):
#             if chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] == pos_arr2[idx2]:
#                 print('NOW')
#                 result = np.abs(meth_arr1[idx1] - meth_arr2[idx2])
#                 idx1 += 1
#                 idx2 += 1
#                 count += 1
#             elif chr_arr1[idx1] < chr_arr2[idx2] or (chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] < pos_arr2[idx2]):
#                 idx1 += 1
#             elif chr_arr1[idx1] > chr_arr2[idx2] or (chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] > pos_arr2[idx2]):
#                 idx2 += 1
#
#     if count > 0:
#         result /= count
#
#     return result
#
#
# def numba_get_pairwise_dissimilarity(args):
#     row_label, col_label = args
#     df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     row_label, col_label = row_label.split('/')[-1].replace('.bed', ''), col_label.split('/')[-1].replace('.bed', '')
#     pair_dist = numba_pair_dist(df1['chr'].values, df2['chr'].values, df1['pos'].values,
#                                 df2['pos'].values, df1['frac_methyl'].values, df2['frac_methyl'].values)
#     print(pair_dist)
#     return row_label, col_label, pair_dist
#
# def numba_pdclust(path_to_bed_files, plot=True):
#     if not path_to_bed_files.endswith('/'):
#         path_to_bed_files = path_to_bed_files + '/'
#
#     file_list = [path_to_bed_files +
#                  f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
#     pool = mp.Pool(processes=mp.cpu_count())
#     # pool = mp.Pool(processes=4)
#
#     results = pool.map(numba_get_pairwise_dissimilarity,
#                        [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])
#
#     index = columns = [x.split('/')[-1].replace('.bed', '') for x in file_list]
#     result_df = pd.DataFrame(index=index, columns=columns)
#
#     for i, elem in enumerate(results):
#         elem1, elem2, pairwise_dissimilarity = elem
#         result_df.at[elem1, elem2] = pairwise_dissimilarity
#     result_df.fillna(0.0, inplace=True)
#     if plot:
#         plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap.png')
#     pool.close()
#     return result_df


import numpy as np
import pandas as pd
import numba as nb
import multiprocessing as mp
import os


# @numba.jit(forceobj=True)
# # numba.f4(numba.i8[:], numba.i8[:], numba.i8[:], numba.i8[:], numba.i8[:], numba.i8[:])
# def numba_pair_dist(chr_arr1, chr_arr2, pos_arr1, pos_arr2, meth_arr1, meth_arr2):
#     common_elements = np.intersect1d(pos_arr1, pos_arr2)
#     print(common_elements)
#     result = 0
#     count = 0
#     idx1 = 0
#     idx2 = 0
#
#     while idx1 < len(chr_arr1) and idx2 < len(chr_arr2):
#         if chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] == pos_arr2[idx2]:
#             result = result + np.abs(meth_arr1[idx1] - meth_arr2[idx2])
#             idx1 = idx1 + 1
#             idx2 = idx2 + 1
#             count = count + 1
#         elif chr_arr1[idx1] < chr_arr2[idx2] or (
#                 chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] < pos_arr2[idx2]):
#             idx1 = idx1 + 1
#         elif chr_arr1[idx1] > chr_arr2[idx2] or (
#                 chr_arr1[idx1] == chr_arr2[idx2] and pos_arr1[idx1] > pos_arr2[idx2]):
#             idx2 = idx2 + 1
#
#     if count > 0:
#         result = result / count
#
#     return result
#
#
# def numba_get_pairwise_dissimilarity(args):
#     row_label, col_label = args
#     df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     row_label, col_label = row_label.split('/')[-1].replace('.bed', ''), col_label.split('/')[-1].replace('.bed', '')
#     df1['chr'] = [x[-1] for x in df1['chr'].values]
#     df2['chr'] = [x[-1] for x in df2['chr'].values]
#     pair_dist = numba_pair_dist(df1['chr'].values, df2['chr'].values, df1['pos'].values,
#                                 df2['pos'].values, df1['frac_methyl'].values, df2['frac_methyl'].values)
#     return row_label, col_label, pair_dist
#
#
# def numba_pdclust(path_to_bed_files, plot=True):
#     if not path_to_bed_files.endswith('/'):
#         path_to_bed_files = path_to_bed_files + '/'
#
#     file_list = [path_to_bed_files + f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
#     pool = mp.Pool(processes=mp.cpu_count())
#
#     results = pool.map(numba_get_pairwise_dissimilarity,
#                        [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])
#
#     index = columns = [x.split('/')[-1].replace('.bed', '') for x in file_list]
#     pairwise_dissimilarity = np.zeros((len(index), len(columns)))
#
#     for i, elem in enumerate(results):
#         elem1, elem2, dissimilarity = elem
#         pairwise_dissimilarity[i, i] = dissimilarity
#
#     result_df = pd.DataFrame(pairwise_dissimilarity, index=index, columns=columns)
#     result_df.fillna(0.0, inplace=True)
#
#     if plot:
#         plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap_numba.png')
#
#     pool.close()
#     return result_df

from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing as mp
from numba import jit

import numpy as np
from numba import njit


import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from numba import njit

@numba.njit(nb.types.Tuple((nb.int64, nb.int64))(nb.int64[:], nb.int64[:], nb.float64[:], nb.float64[:]))
def numba_calc_pdist(pos_arr1, pos_arr2, meth_arr1, meth_arr2):
    idx1 = 0
    idx2 = 0
    number_of_cpgs = 0
    pair_dist = 0
    while idx1 < len(pos_arr1) and idx2 < len(pos_arr2):
        if pos_arr1[idx1] == pos_arr2[idx2]:
            number_of_cpgs = number_of_cpgs + 1
            pair_dist = pair_dist + np.abs(meth_arr1[idx1] - meth_arr2[idx2])
        elif pos_arr1[idx1] < pos_arr2[idx2]:
            idx1 = idx1 + 1
        else:
            idx2 = idx2 + 1
    return pair_dist, number_of_cpgs

def numba_get_pairwise_dissimilarity(args):
    total_cpgs = 0
    pairwise_dissimilarity = 0
    row_label, col_label = args
    df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    row_label, col_label = Path(row_label).stem, Path(col_label).stem
    for i in range(0, 20):
        tmp_df1 = df1[df1['chr'] == f'chr{i}']
        tmp_df2 = df2[df2['chr'] == f'chr{i}']
        pair_dist, cpgs = numba_calc_pdist(tmp_df1['pos'].to_numpy(), tmp_df2['pos'].to_numpy(), tmp_df1['frac_methyl'].to_numpy(), tmp_df2['frac_methyl'].to_numpy())
        pairwise_dissimilarity = pairwise_dissimilarity + pair_dist
        total_cpgs = total_cpgs + cpgs
    return row_label, col_label, (pairwise_dissimilarity / total_cpgs)


def numba_pdclust(path_to_bed_files, plot=True):
    path_to_bed_files = Path(path_to_bed_files)
    file_list = list(path_to_bed_files.glob("*.bed"))

    with mp.Pool(processes=4) as pool:
        results = pool.map(numba_get_pairwise_dissimilarity, [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for elem in results:
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity

    result_df.fillna(0.0, inplace=True)

    if plot:
        plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap_numba.png')

    return result_df


######################################
########### withouth numba ###########
######################################

# def get_pairwise_dissimilarity(args):
#     """
#     Calculate pairwise dissimilarity between two cells :param args: tuple of two .bed files where each bed file
#     contains the columns 'chromosome', 'start', 'end', 'methylation' (in percent, e.g. 0 or 100) :return: tuple of (
#     row_label, col_label, pairwise_dissimilarity) where row_label is the name of the first cell, col_label is the
#     name of the second cell, and pairwise_dissimilarity is the hamming distance between the two cells
#     """
#
#     row_label, col_label = args
#     df1 = pd.read_csv(row_label, sep='\t', usecols=(
#         0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     df2 = pd.read_csv(col_label, sep='\t', usecols=(
#         0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
#     row_label, col_label = row_label.split(
#         '/')[-1].replace('.bed', ''), col_label.split('/')[-1].replace('.bed', '')
#     merged_df = pd.merge(df1, df2, on=['chr', 'pos'], suffixes=(
#         '_df1', '_df2'), how='inner')
#     total_cpgs = merged_df.shape[0]
#     pairwise_dissimilarity_total = np.sum(
#         np.abs(merged_df['frac_methyl_df1'] - merged_df['frac_methyl_df2']))
#     if total_cpgs != 0:
#         pairwise_dissimilarity = pairwise_dissimilarity_total / total_cpgs
#     else:
#         pairwise_dissimilarity = 0.0
#     return row_label, col_label, pairwise_dissimilarity
#
#
# def trivial_pdclust(path_to_bed_files, plot=True):
#     if not path_to_bed_files.endswith('/'):
#         path_to_bed_files = path_to_bed_files + '/'
#
#     file_list = [path_to_bed_files +
#                  f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
#     # pool = mp.Pool(processes=mp.cpu_count())
#     pool = mp.Pool(processes=4)
#
#     results = pool.map(get_pairwise_dissimilarity,
#                        [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])
#
#     index = columns = [x.split('/')[-1].replace('.bed', '') for x in file_list]
#     result_df = pd.DataFrame(index=index, columns=columns)
#
#     for i, elem in enumerate(results):
#         elem1, elem2, pairwise_dissimilarity = elem
#         result_df.at[elem1, elem2] = pairwise_dissimilarity
#     result_df.fillna(0.0, inplace=True)
#     if plot:
#         plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap_trivial.png')
#     pool.close()
#     return result_df


from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing as mp


def get_pairwise_dissimilarity(args):
    row_label, col_label = args
    df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    row_label, col_label = Path(row_label).stem, Path(col_label).stem
    merged_df = pd.merge(df1, df2, on=['chr', 'pos'], suffixes=('_df1', '_df2'), how='inner')
    total_cpgs = merged_df.shape[0]
    pairwise_dissimilarity_total = np.sum(np.abs(merged_df['frac_methyl_df1'] - merged_df['frac_methyl_df2']))
    pairwise_dissimilarity = pairwise_dissimilarity_total / total_cpgs if total_cpgs != 0 else 0.0
    return row_label, col_label, pairwise_dissimilarity


def trivial_pdclust(path_to_bed_files, plot=True):
    path_to_bed_files = Path(path_to_bed_files)
    file_list = list(path_to_bed_files.glob("*.bed"))

    with mp.Pool(processes=4) as pool:
        results = pool.map(get_pairwise_dissimilarity, [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for elem in results:
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity

    result_df.fillna(0.0, inplace=True)

    if plot:
        plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap_trivial.png')

    return result_df


##########################################
############# with dask ##################
##########################################

import dask.dataframe as dd
import dask.distributed
from dask.distributed import wait

def dask_get_pairwise_dissimilarity(args):
    row_label, col_label = args
    df1 = dd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'],
                      dtype={'pos': 'int64', 'frac_methyl': 'float64'})
    df2 = dd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'],
                      dtype={'pos': 'int64', 'frac_methyl': 'float64'})
    row_label, col_label = row_label.split('/')[-1].replace('.bed', ''), col_label.split('/')[-1].replace('.bed', '')
    merged_df = dd.merge(df1, df2, on=['chr', 'pos'], suffixes=('_df1', '_df2'), how='inner')
    total_cpgs = merged_df.shape[0]
    pairwise_dissimilarity_total = merged_df['frac_methyl_df1'].sub(merged_df['frac_methyl_df2']).abs().sum()
    pairwise_dissimilarity = pairwise_dissimilarity_total / total_cpgs
    return row_label, col_label, float(pairwise_dissimilarity.compute())


def dask_pdclust(path_to_bed_files, plot=True):
    if not path_to_bed_files.endswith('/'):
        path_to_bed_files = path_to_bed_files + '/'

    file_list = [path_to_bed_files + f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]

    cluster = dask.distributed.LocalCluster()  # Erstellt einen lokalen Cluster
    client = dask.distributed.Client(cluster)  # Verbindet den Client mit dem Cluster

    results = []
    for elem1 in file_list:
        for elem2 in file_list:
            if elem1 != elem2:
                results.append(client.submit(dask_get_pairwise_dissimilarity, (elem1, elem2)))

    wait(results)  # Warten Sie auf das Abschließen aller Berechnungen

    index = columns = [x.split('/')[-1].replace('.bed', '') for x in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for i, elem in enumerate(results):
        elem1, elem2, pairwise_dissimilarity = elem.result()
        result_df.at[elem1, elem2] = pairwise_dissimilarity

    result_df.fillna(0.0, inplace=True)

    if plot:
        plot_heatmap(result_df, output_path=f'{path_to_bed_files}/heatmap_dask.png')

    client.close()  # Schließt die Verbindung zum Cluster
    cluster.close()  # Beendet den lokalen Cluster

    return result_df




##########################################


def plot_heatmap(dist_matrix, output_path='./heatmap.png'):
    """
    Plot heatmap of pairwise dissimilarity between all cells with annotation bars depending on cell types
    :param output_path: path to save the heatmap
    :param dist_matrix: return value from pdclust(). Alternatively: Distance matrix between all cells (pandas dataframe)
    :return: None, saves heatmap at the destination of output_path
    """
    # colors = ['red' if 'ESLAM' in x else 'blue' for x in dist_matrix.columns]
    cluster_map = sns.clustermap(data=dist_matrix, method='ward', metric='euclidean', figsize=(20, 20), cmap='YlOrRd',
                                 xticklabels=False, yticklabels=False, vmin=8.5, vmax=12.5)

    cluster_map.ax_row_dendrogram.set_visible(False)
    cluster_map.ax_col_dendrogram.set_visible(True)

    plt.savefig(output_path)


import timeit

if __name__ == '__main__':
    # df_trivial = trivial_pdclust('/home/ptheo/Documents/university/bachelor/data/human_data/test/')
    # print(df_trivial)
    # df_numba = numba_pdclust('/home/ptheo/Documents/university/bachelor/data/human_data/test/')
    # print(df_numba)
    # print('Trivial:' + str(timeit.timeit("trivial_pdclust", globals=locals(), number=2000000)))
    # print('Numba:' + str(timeit.timeit("numba_pdclust", globals=locals(), number=2000000)))
    # print('Dask:' + str(timeit.timeit("dask_pdclust", globals=locals(), number=2000000)))
    df_dask = dask_pdclust('/home/ptheo/Documents/university/bachelor/data/human_data/')
    df_dask.to_csv('df_dask.csv')


def make_unique_combinations(path_to_bed_files):
    """
    Make unique combinations of all cells in a directory
    :param path_to_bed_files: containing the BED files
    :return: list of tuples of unique combinations of all cells
    """
    if not path_to_bed_files.endswith('/'):
        path_to_bed_files = path_to_bed_files + '/'
    file_list = [path_to_bed_files +
                 f for f in os.listdir(path_to_bed_files) if f.endswith('.bed')]
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
