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

def get_nonempty_accessions(metadata, g4_list, nonempty_result):
    with open(g4_list, 'r') as file:
        file_list = file.readlines()
        file_list = [file.strip() for file in file_list]
                
    # read metadata file
    df_meta = pd.read_csv(metadata, sep="\t")
    df_meta['filename'] = df_meta['file'].apply(lambda x: x.split("/")[-1])
    df_meta['type'] = df_meta['file'].apply(lambda x: get_type(x))
    df_meta = df_meta[~df_meta['type'].isin(["cds", "rna"])]
    file_size_dict = dict(zip(df_meta['filename'], df_meta['empty']))
    
    # keep only files with data
    df = pd.DataFrame(file_list, columns=["filename"])
    df["has_data"] = df["filename"].apply(lambda x: file_size_dict[x] if x in file_size_dict else True)
    df = df[df["has_data"]]
    df.drop(["has_data"], axis=1, inplace=True)
    print("Non-empty files:", len(df))
    df.to_csv(nonempty_result, sep="\t", index=False, header=False)

if __name__ == '__main__':
    # Get list of unique accessions from g4hunter
    g4hunter_metadata = f"{BASE_PATH}/metadata/g4_file_metadata.txt"
    g4hunter_list = f"{BASE_PATH}/slurm/files/g4_list.txt"
    g4hunter_nonempty_list = f"{BASE_PATH}/slurm/files/g4_nonempty_list.txt"
    get_nonempty_accessions(g4hunter_metadata, g4hunter_list, g4hunter_nonempty_list)

    # Get list of unique accessions from regex
    regex_metadata = f"{BASE_PATH}/metadata/regex_file_metadata.txt"
    regex_list = f"{BASE_PATH}/slurm/files/regex_list.txt"
    regex_nonempty_list = f"{BASE_PATH}/slurm/files/regex_nonempty_list.txt"
    get_nonempty_accessions(regex_metadata, regex_list, regex_nonempty_list)