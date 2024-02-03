import pandas as pd

BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

g4hunter_path = f"{BASE_PATH}/raw_data/g4hunter"
g4hunter_nonempty_list = f"{BASE_PATH}/slurm/files/g4_nonempty_list.txt"
g4hunter_result_path = f"{BASE_PATH}/data/g4hunter_bed"
df_g4hunter = pd.read_csv(g4hunter_nonempty_list, sep="\t", header=None)
df_g4hunter.columns = ["filename"]
df_g4hunter["path"] = g4hunter_path
df_g4hunter["result_path"] = g4hunter_result_path

regex_path = f"/storage/group/izg5139/default/external/genomes/all_genomes/regex/regex_gresults"
regex_nonempty_list = f"{BASE_PATH}/slurm/files/regex_nonempty_list.txt"
regex_result_path = f"{BASE_PATH}/data/regex_bed"
df_regex = pd.read_csv(regex_nonempty_list, sep="\t", header=None)
df_regex.columns = ["filename"]
df_regex["path"] = regex_path
df_regex["result_path"] = regex_result_path

df = pd.concat([df_g4hunter, df_regex])
# reorder columns to match the order in which they are passed to the script
df = df[["path", "filename", "result_path"]]
# save as txt file with space as delimiter
df.to_csv(f"{BASE_PATH}/slurm/files/csv_to_bed_args.txt", sep=" ", index=False, header=False)