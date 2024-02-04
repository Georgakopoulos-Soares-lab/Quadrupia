import os
import json
import pandas as pd

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

# combine results from g4_distribution 
g4_hunter_results = []
regex_results = []
files = os.listdir(f"{BASE_PATH}/data/g4_distribution")
l = len(files)
for i in range(l):
    file = files[i]
    with open(f"{BASE_PATH}/data/g4_distribution/{file}") as f:
        data = json.load(f)
    g4_hunter_results.append(data['g4hunter'])
    regex_results.append(data['regex'])

# save results
df_g4hunter = pd.DataFrame(g4_hunter_results)
df_g4hunter.to_csv(f"{BASE_PATH}/results/g4hunter_g4_distribution.csv", index=False)
df_regex = pd.DataFrame(regex_results)
df_regex.to_csv(f"{BASE_PATH}/results/regex_g4_distribution.csv", index=False)