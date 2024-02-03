import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"
# Path to directory containing GFF files
GFF_DIR_PATH = "/storage/group/izg5139/default/external/gff"

if __name__ == '__main__':
    regex_list = f"{BASE_PATH}/slurm/files/regex_list.txt"
    
    # Only keep items in gff_files that are in regex_list
    gff_files = os.listdir(GFF_DIR_PATH)
    gff_files = [file.replace('gff', 'csv') for file in gff_files]
    g4s_list = open(regex_list).readlines()
    g4s_list = [file.strip() for file in g4s_list]
    print(len(gff_files))
    gff_files = [file.replace('csv', 'gff') for file in gff_files if file in g4s_list]
    print(len(gff_files))
    with open(f"{BASE_PATH}/slurm/files/gff_files.txt", 'w') as f:
        for file in gff_files:
            f.write(file + '\n')