import glob
import shutil
import pandas as pd

files = glob.glob("/home/ptheo/Documents/university/bachelor/data/human_data/*_filtered.bed")
chromosomes = ['chr' + str(i) for i in range(1, 20)]

# if __name__ == '__main__':
#     for file in files:
#         name = file.split('/')[-1].split('_')[-1].replace('.cpg.bed', '')
# 
#         # rename files to remove unneccessary information in file name
#         shutil.copyfile(file, f"/home/ptheo/Documents/university/bachelor/data/human_data/{name}.bed")
# 
#         # remove rows from bed file to keep only information of autosomal cells
#         bed = pd.read_csv(file, sep='\t', names=['chr', 'start', 'end', 'frac_methyl', 'a', 'b', 'c'])
#         bed = bed.drop(bed[bed['chr'] not in chromosomes])
#         bed = bed[bed['chr'].isin(chromosomes)]
#         bed.to_csv(f'/home/ptheo/Documents/university/bachelor/data/human_data/{name}_filtered.bed', sep='\t')
# 
#         # delete all cells with CpG count < 130.000


if __name__ == '__main__':
    delete = []
    bed_files = glob.glob('/home/ptheo/Documents/university/bachelor/data/human_data/filtered/*.bed')
    for file in bed_files:
        with open(file, "rbU") as f:
            num_lines = sum(1 for _ in f)
        if num_lines < 130000:
            delete.append(file)

    print(delete)