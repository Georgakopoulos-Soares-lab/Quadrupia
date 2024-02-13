import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

G4_PATH = f"{BASE_PATH}/raw_data/g4hunter"
G4_LIST = f"{BASE_PATH}/slurm/files/g4_list.txt"

REGEX_PATH = f"{BASE_PATH}/raw_data/regex"
REGEX_LIST = f"{BASE_PATH}/slurm/files/regex_list.txt"
    
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
    ## get a list of unique accessions from G4Hunter
    
    # Get list of accessions in G4Hunter
    file_list = os.listdir(G4_PATH)
    df = pd.DataFrame(file_list, columns=["filename"])
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
    
    ## get a list of regex accessions that are not in g4hunter
    
    # remove regex accessions already found in g4hunter
    g4_list = df["filename"].tolist()
    g4_list_alt = [('GCF' + file[3:]) for file in g4_list] + [('GCA' + file[3:]) for file in g4_list]
    g4_combined_list = g4_list + g4_list_alt
    regex_files = os.listdir(REGEX_PATH)
    regex_extra = list(set(regex_files) - set(g4_combined_list))
    print("Files in regex but not in G4Hunter:", len(regex_extra))
    
    # remove duplicate accessions and keep only GCF accession if both GCF and GCA are present
    df_reg_extra = pd.DataFrame(regex_extra, columns=["filename"])
    df_reg_extra["accession"] = df_reg_extra["filename"].apply(lambda x: '_'.join(x.split("_")[1:]))
    df_reg_extra = df_reg_extra.groupby("accession").agg({
        "filename": lambda x: prefer_GCF(x)
    }).reset_index()
    regex_extra_list = df_reg_extra["filename"].tolist()
    print("Extra regex files after removing duplicate accessions:", len(regex_extra_list))
    print(regex_extra_list)
    
    regex_final_list = regex_extra_list + g4_list
    with open(REGEX_LIST, 'w') as f:
        for file in regex_final_list:
            f.write(file + '\n')