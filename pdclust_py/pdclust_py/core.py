from pathlib import Path
import numba
import pandas as pd
import numpy as np
import multiprocessing as mp
from tqdm import tqdm

@numba.jit(nopython=True)
def numba_calc_pdist(pos_arr1, pos_arr2, meth_arr1, meth_arr2):
    """
    Calculate the pairwise distance between two cells
    :param pos_arr1: Array of starting positions of CpGs in the first cell
    :param pos_arr2: Array of starting positions of CpGs in the second cell
    :param meth_arr1: Array of methylation values in the first cell
    :param meth_arr2: Array of methylation values in the second cell
    :return: pair_dist = pairwise distance, number_of_cpgs = Number of CpGs that are present in both cells
    """
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
    """
    Open the bed files and pass the data to numba_calc_pdist to calculate the pairwise distance
    :param args: Tuple of two bed file paths
    :return: row_label = name of cell 1, col_label = name of cell 2, pairwise_dissimilarity
    """
    total_cpgs = 0
    pairwise_dissimilarity = 0
    row_label, col_label = args
    df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'])
    df1 = df1.astype({'pos': np.int32, 'frac_methyl': np.int8})
    df2 = df2.astype({'pos': np.int32, 'frac_methyl': np.int8})
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


def numba_pdclust(path_to_bed_files):
    """
    Calculate the pairwise dissimilarity between all cells in the given directory
    :param path_to_bed_files: Directory containing the bed files
    :return: Pandas DataFrame containing all pairwise dissimilarities in square form
    """
    path_to_bed_files = Path(path_to_bed_files)
    file_list = list(path_to_bed_files.glob("*.bed"))
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(numba_get_pairwise_dissimilarity,
                           [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2])

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)
    for elem in results:
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity
    result_df = result_df[result_df.index.to_list()]
    result_df = result_df.fillna(0.0)
    return result_df


#########################################
############# end numba #################
#########################################

def get_pairwise_dissimilarity(tuple_of_bed_file_paths):
    """
    Open the bed files and calculate the pairwise distance
    :param tuple_of_bed_file_paths: Tuple of two bed file paths
    :return: row_label = name of cell 1, col_label = name of cell 2, pairwise_dissimilarity
    """
    row_label, col_label = tuple_of_bed_file_paths
    df1 = pd.read_csv(row_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'], dtype={'chr': 'string', 'pos': np.int16, 'frac_methyl': np.int8})
    df2 = pd.read_csv(col_label, sep='\t', usecols=(0, 1, 3), names=['chr', 'pos', 'frac_methyl'], dtype={'chr': 'string', 'pos': np.int16, 'frac_methyl': np.int8})
    row_label, col_label = Path(row_label).stem, Path(col_label).stem
    merged_df = pd.merge(df1, df2, on=['chr', 'pos'], suffixes=('_df1', '_df2'), how='inner')
    total_cpgs = merged_df.shape[0]
    pairwise_dissimilarity_total = np.sum(np.abs(merged_df['frac_methyl_df1'] - merged_df['frac_methyl_df2']))
    pairwise_dissimilarity = pairwise_dissimilarity_total / total_cpgs if total_cpgs != 0 else 0.0
    return row_label, col_label, pairwise_dissimilarity


def trivial_pdclust(path_to_bed_files):
    """
    :param path_to_bed_files: Path to directory containing the bed files
    :return: Pandas DataFrame containing all pairwise dissimilarities in square form
    """
    path_to_bed_files = Path(path_to_bed_files)
    file_list = sorted(list(path_to_bed_files.glob("*.bed")))

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(get_pairwise_dissimilarity, [(elem1, elem2) for elem1 in file_list for elem2 in file_list if elem1 != elem2]), total=len(file_list) * (len(file_list) - 1) / 2))

    index = columns = [file.stem for file in file_list]
    result_df = pd.DataFrame(index=index, columns=columns)
    for elem in results:
        elem1, elem2, pairwise_dissimilarity = elem
        result_df.at[elem1, elem2] = pairwise_dissimilarity
    result_df = result_df[result_df.index.to_list()]
    result_df = result_df.fillna(0.0)
    
    return result_df

