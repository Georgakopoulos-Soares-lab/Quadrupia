import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"
# Path to directory containing GFF files
GFF_DIR_PATH = "/storage/group/izg5139/default/external/gff"

def get_accession(file):
    return '_'.join(file.strip().split('_')[:2])

if __name__ == '__main__':
    g4_list = f"{BASE_PATH}/slurm/files/g4_list.txt"
    
    # Only keep items in gff_files that are in g4_list
    gff_files = os.listdir(GFF_DIR_PATH)
    g4s_list = open(g4_list).readlines()
    g4s_list = [get_accession(file) for file in g4s_list]
    print(len(gff_files))
    # get GFF files that have a corresponding G4 file
    gff_files = [file for file in gff_files if get_accession(file) in g4s_list]
    print(len(gff_files))
    with open(f"{BASE_PATH}/slurm/files/gff_files.txt", 'w') as f:
        for file in gff_files:
            f.write(file + '\n')