import os
import json
import traceback
import pandas as pd
from multiprocessing import Process, Lock, Manager

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

# GFF files(in BED format) for getting gene data
GFF_BED_PATH = f"{BASE_PATH}/data/gff_bed"

# Path to result file
RESULT_PATH = f"{BASE_PATH}/results/gff_feature_data.csv"
      
# number of CPUs on the machine
NUM_PROC = os.cpu_count()
print(NUM_PROC)
    
def get_gff_data(filename, gff_feature_data, lock):
    print(f"Getting gff data for {filename}")
    try:
        gff_bed_file = os.path.join(GFF_BED_PATH, filename)
        gff = pd.read_csv(gff_bed_file, sep="\t", header=None)
        gff.columns = ['chr', 'start', 'end', 'type']
        gff['length'] = gff['end'] - gff['start']
        gff_features = gff.groupby('type').agg({'length': 'sum'}).reset_index()
        # create a dictionary with feature type as key and length as value
        row = gff_features.set_index('type').T.to_dict('records')[0]
        row['accession'] = '_'.join(filename.split('_')[:2])
        lock.acquire()
        gff_feature_data.append(row)
        lock.release()
        
    except Exception as e:
        print(traceback.format_exc())
        
if __name__ == "__main__":
    # get list of gff bed files
    gff_file_list = os.listdir(GFF_BED_PATH)
    
    # create temp directory if not exists
    if not os.path.exists(f"{BASE_PATH}/temp"):
        os.makedirs(f"{BASE_PATH}/temp")   
    
    # shared variable between processes
    gff_feature_data = Manager().list()
    
    # create a lock for the shared variable
    lock = Lock()

    # create processes for each G4 file
    index = 0
    while index < len(gff_file_list):
        processes = []
        for _ in range(NUM_PROC):
            if index == len(gff_file_list):
                break
            process = Process(target=get_gff_data, args=(gff_file_list[index], gff_feature_data, lock))
            processes.append(process)
            index += 1
            
        # start processes
        for process in processes:
            process.start()

        # wait for all processes to finish
        # block the main programm until these processes are finished
        for process in processes:
            process.join()
        
    # save the results to file after all processes are executed 
    df_result = pd.DataFrame(list(gff_feature_data))
    df_result.to_csv(RESULT_PATH, index=False)
