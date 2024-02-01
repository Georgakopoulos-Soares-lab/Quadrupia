import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

G4_PATH = f"{BASE_PATH}/raw_data/g4hunter"
G4_LIST = f"{BASE_PATH}/slurm/files/g4_list.txt"

REGEX_PATH = f"/storage/group/izg5139/default/external/genomes/all_genomes/regex/regex_gresults"
REGEX_LIST = f"{BASE_PATH}/slurm/files/regex_list.txt"

def get_type(filename):
    if "cds_from_genomic" in filename:
        return "cds"
    elif "rna_from_genomic" in filename:
        return "rna"
    else:
        return "genomic"
    
def prefer_GCF(files):
    files = list(files)
    if len(files) == 1:
        return files[0]
    else:
        for file in files:
            if file[:3] == "GCF":
                return file
        return files[0]

if __name__ == '__main__':
    # Get list of unique accessions in G4Hunter
    file_list = os.listdir(G4_PATH)
    
    # process only genomic file
    files = [file for file in file_list \
        if ("_genomic" in file) \
            and ("cds_from_genomic" not in file) \
                and ("rna_from_genomic" not in file)]
    
    # read file list
    df = pd.DataFrame(files, columns=["filename"])
    df["accession"] = df["filename"].apply(lambda x: '_'.join(x.split("_")[1:]))
    print("All files:", len(df))
    
    # keep only GCF accession if both GCF and GCA are present
    # sort by accession number and accession in descending order and keep only the first row (GCF) 
    # group by accession number
    df = df.groupby("accession").agg({
        "filename": lambda x: prefer_GCF(x)
    }).reset_index()
    df.drop(["accession"], axis=1, inplace=True)
    df.sort_values("filename", ascending=True, inplace=True)
    
    print("G4Hunter files after removing duplicate accessions:", len(df))
    df.to_csv(G4_LIST, sep="\t", index=False, header=False)
    
    # find accessions in regex that are not in g4hunter
    g4_list = df["filename"].tolist()
    regex_files = os.listdir(REGEX_PATH)
    regex_files = [file for file in regex_files \
        if ("_genomic" in file) \
            and ("cds_from_genomic" not in file) \
                and ("rna_from_genomic" not in file)]
    
    regex_extra = [file for file in regex_files \
        if file.replace('GCA', 'GCF') not in g4_list and file.replace('GCF', 'GCA') not in g4_list]
    
    df_reg_extra = pd.DataFrame(regex_extra, columns=["filename"])
    df_reg_extra["accession"] = df_reg_extra["filename"].apply(lambda x: '_'.join(x.split("_")[1:]))
    print("Files in regex but not in G4Hunter:", len(df_reg_extra))
    
    # remove duplicate accessions and keep only GCF accession if both GCF and GCA are present
    df_reg_extra = df_reg_extra.groupby("accession").agg({
        "filename": lambda x: prefer_GCF(x)
    }).reset_index()
    df_reg_extra.drop(["accession"], axis=1, inplace=True)
    regex_list = df_reg_extra["filename"].tolist()
    print("Extra regex files after removing duplicate accessions:", len(regex_list))
    print(regex_list)
    
    regex_final_list = regex_list + g4_list
    with open(REGEX_LIST, 'w') as f:
        for file in regex_final_list:
            f.write(file + '\n')