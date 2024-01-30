import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"
# Path to directory containing GFF files
GFF_DIR_PATH = "/storage/group/izg5139/default/external/gff"

if __name__ == '__main__':
    
    # Create a list which is a union of g4hunter and regex file lists
    g4hunter_list = f"{BASE_PATH}/slurm/files/g4hunter_list.txt"
    regex_list = f"{BASE_PATH}/slurm/files/regex_list.txt"
    combined_list = f"{BASE_PATH}/slurm/files/combined_list.txt"
    cmd = f"cat {g4hunter_list} {regex_list} | sort | uniq > {combined_list}"
    os.system(cmd)
    
    # Only keep items in gff_files that are in combined_list
    gff_files = os.listdir(GFF_DIR_PATH)
    gff_files = [file.replace('gff', 'csv') for file in gff_files]
    combined_g4s_list = open(combined_list).readlines()
    combined_g4s_list = [file.strip() for file in combined_g4s_list]
    print(len(gff_files))
    gff_files = [file for file in gff_files if file in combined_g4s_list]
    print(len(gff_files))
    with open(f"{BASE_PATH}/slurm/files/gff_files.txt", 'w') as f:
        for file in gff_files:
            f.write(file + '\n')