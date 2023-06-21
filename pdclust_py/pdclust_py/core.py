import os
import itertools
import numba
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import multiprocessing as mp

from pathlib import Path
from tqdm import tqdm


@numba.jit(nopython=True)
def numba_calc_pdist(pos_arr1, pos_arr2, meth_arr1, meth_arr2):
    idx1 = 0
    idx2 = 0
    number_of_cpgs = 0
    pair_dist = 0
    while idx1 < len(pos_arr1) and idx2 < len(pos_arr2):
        if pos_arr1[idx1] == pos_arr2[idx2]:
            number_of_cpgs = number_of_cpgs + 1
            pair_dist = pair_dist + np.abs(meth_arr1[idx1] - meth_arr2[idx2])
            idx1 = idx1 + 1
            idx2 = idx2 + 1
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
    df1 = df1.astype({'pos': 'i8', 'frac_methyl': 'i4'})
    row_label, col_label = Path(row_label).stem, Path(col_label).stem
    for i in range(0, 20):
        tmp_df1 = df1[df1['chr'] == f'chr{i}'].copy(deep=True)
        tmp_df2 = df2[df2['chr'] == f'chr{i}'].copy(deep=True)
        pair_dist, cpgs = numba_calc_pdist(tmp_df1['pos'].to_numpy(), tmp_df2['pos'].to_numpy(),
                                           tmp_df1['frac_methyl'].to_numpy(),
                                           tmp_df2['frac_methyl'].to_numpy())
        pairwise_dissimilarity = pairwise_dissimilarity + pair_dist
        total_cpgs = total_cpgs + cpgs
    return row_label, col_label, (pairwise_dissimilarity / total_cpgs)


def numba_pdclust(path_to_bed_files, plot=True):
    results = []
    path_to_bed_files = Path(path_to_bed_files)
    file_list = list(path_to_bed_files.glob("*.bed"))
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(numba_get_pairwise_dissimilarity,
                           [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for i, elem in enumerate(results):
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity
        
    result_df.fillna(0.0, inplace=True)
    
    return result_df


#########################################
############# end numba #################
#########################################

import numpy as np
import os
import tqdm
import seaborn as sns
import datetime
import multiprocessing as mp


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

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(get_pairwise_dissimilarity, [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)

    for elem in results:
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity

    result_df.fillna(0.0, inplace=True)

    return result_df
