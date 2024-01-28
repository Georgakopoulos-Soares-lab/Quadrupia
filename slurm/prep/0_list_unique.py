import os
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

def get_type(filename):
    if "cds_from_genomic" in filename:
        return "cds"
    elif "rna_from_genomic" in filename:
        return "rna"
    else:
        return "genomic"

def get_unique_accessions(filepath, metadata, result, nonempty_result):
    file_list = os.listdir(filepath)
    
    # process only genomic file
    files = [file for file in file_list \
        if ("_genomic" in file) \
            and ("cds_from_genomic" not in file) \
                and ("rna_from_genomic" not in file)]
                
    # read metadata file
    df_meta = pd.read_csv(metadata, sep="\t")
    df_meta['filename'] = df_meta['file'].apply(lambda x: x.split("/")[-1])
    df_meta['type'] = df_meta['file'].apply(lambda x: get_type(x))
    df_meta = df_meta[~df_meta['type'].isin(["cds", "rna"])]
    file_size_dict = dict(zip(df_meta['filename'], df_meta['empty']))
    
    # read file list
    df = pd.DataFrame(files, columns=["filename"])
    df["accession"] = df["filename"].apply(lambda x: '_'.join(x.split("_")[:2]))
    df["accession_number"] = df["accession"].apply(lambda x: x.split("_")[1])
    print("All files:", len(df))
    
    # keep only GCF accession if both GCF and GCA are present
    # sort by accession number and accession in descending order and keep only the first row (GCF) 
    df.sort_values(["accession_number", "accession"], ascending=False, inplace=True)
    df.drop_duplicates(subset=["accession_number"], keep="first", inplace=True)
    df.drop(["accession", "accession_number"], axis=1, inplace=True)
    print("Files after removing duplicate accessions:", len(df))
    df.sort_values("filename", ascending=True, inplace=True)
    df.to_csv(result, sep="\t", index=False, header=False)
    
    # keep only files with data
    df["has_data"] = df["filename"].apply(lambda x: file_size_dict[x] if x in file_size_dict else True)
    df = df[df["has_data"]]
    df.drop(["has_data"], axis=1, inplace=True)
    print("Non-empty files:", len(df))
    df.to_csv(nonempty_result, sep="\t", index=False, header=False)

if __name__ == '__main__':
    # Get list of unique accessions from g4hunter
    g4hunter_path = f"{BASE_PATH}/raw_data/g4hunter"
    g4hunter_metadata = f"{BASE_PATH}/metadata/g4_file_metadata.txt"
    g4hunter_list = f"{BASE_PATH}/slurm/files/g4hunter_list.txt"
    g4hunter_nonempty_list = f"{BASE_PATH}/slurm/files/g4hunter_nonempty_list.txt"
    get_unique_accessions(g4hunter_path, g4hunter_metadata, g4hunter_list, g4hunter_nonempty_list)

    # Get list of unique accessions from regex
    regex_path = f"/storage/group/izg5139/default/external/genomes/all_genomes/regex/regex_gresults"
    regex_metadata = f"{BASE_PATH}/metadata/regex_file_metadata.txt"
    regex_list = f"{BASE_PATH}/slurm/files/regex_list.txt"
    regex_nonempty_list = f"{BASE_PATH}/slurm/files/regex_nonempty_list.txt"
    get_unique_accessions(regex_path, regex_metadata, regex_list, regex_nonempty_list)