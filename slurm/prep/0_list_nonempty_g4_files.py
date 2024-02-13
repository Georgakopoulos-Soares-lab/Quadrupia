import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

def get_nonempty_accessions(empty_list, g4_list, nonempty_result):
    with open(g4_list, 'r') as file:
        file_list = file.readlines()
        file_list = [file.strip() for file in file_list]
    file_set = set(file_list)
                
    # read empty file list
    df_meta = pd.read_csv(empty_list, sep="\t", header=None)
    # keep only 2nd column
    df_meta = df_meta.iloc[:, 1:2]
    df_meta.columns = ["file"]
    df_meta['filename'] = df_meta['file'].apply(lambda x: x.split("/")[-1])
    empty_set = set(df_meta['filename'].tolist())
    
    # get nonempty files
    nonempty_set = file_set - empty_set
    nonempty_list = list(nonempty_set)
    with open(nonempty_result, 'w') as file:
        for line in nonempty_list:
            file.write(f"{line}\n")
    print(f"Nonempty files: {len(nonempty_list)}")
    print(f"Empty files: {len(file_set) - len(nonempty_list)}")

if __name__ == '__main__':
    # Get list of unique accessions from g4hunter
    g4hunter_metadata = f"{BASE_PATH}/metadata/empty_g4_files.txt"
    g4hunter_list = f"{BASE_PATH}/slurm/files/g4_list.txt"
    g4hunter_nonempty_list = f"{BASE_PATH}/slurm/files/g4_nonempty_list.txt"
    get_nonempty_accessions(g4hunter_metadata, g4hunter_list, g4hunter_nonempty_list)

    # Get list of unique accessions from regex
    regex_metadata = f"{BASE_PATH}/metadata/empty_regex_files.txt"
    regex_list = f"{BASE_PATH}/slurm/files/regex_list.txt"
    regex_nonempty_list = f"{BASE_PATH}/slurm/files/regex_nonempty_list.txt"
    get_nonempty_accessions(regex_metadata, regex_list, regex_nonempty_list)