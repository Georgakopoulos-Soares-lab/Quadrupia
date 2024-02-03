import os
import json

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad/data"
OUTPUT_PATH = "/storage/group/izg5139/default/akshatha/gquad/slurm/files"
# Path to directory containing GFF files
GFF_BED_PATH = f"{BASE_PATH}/gff_bed"
    
print("Saving GFF bed file map")
# list of gff bed files
gff_bed_files = os.listdir(GFF_BED_PATH)
# convert list to dictionary with accession as key and gff bed file as value for efficient search
gff_bed_files = {'_'.join(file.split('_')[:2]): file for file in gff_bed_files}
# save dictionary to json file
with open(f"{OUTPUT_PATH}/gff_bed_files.json", 'w') as file:
    json.dump(gff_bed_files, file)