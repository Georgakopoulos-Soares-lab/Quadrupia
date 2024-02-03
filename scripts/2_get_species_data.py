import os
import json
import traceback
import pandas as pd
from multiprocessing import Process, Lock, Manager

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"

# List of genome files
FILE_LIST = f"{BASE_PATH}/slurm/files/regex_list.txt"

# Path to result file
RESULT_PATH = f"{BASE_PATH}/results/species_data.csv"

# Assembly summary files for getting genome details and taxid
GENBANK_SUMMARY = f"{BASE_PATH}/metadata/assembly_summary_genbank.txt"
REFSEQ_SUMMARY = f"{BASE_PATH}/metadata/assembly_summary_refseq.txt"

# Tree of life data for getting taxonomic details
TREE_OF_LIFE = f"{BASE_PATH}/metadata/tree_of_life.csv"

# GFF files(in BED format) for getting gene data
GFF_BED_PATH = "/storage/group/izg5139/default/akshatha/gquad/data/gff_bed"
      
# number of CPUs on the machine
NUM_PROC = os.cpu_count()
print(NUM_PROC)

# read GenBank summary file as csv excluding the first 1 line
df_genbank = pd.read_csv(GENBANK_SUMMARY, sep="\t", skiprows=1)
df_genbank.rename(columns={df_genbank.columns[0]: "assembly_accession"}, inplace=True)
df_genbank.set_index('assembly_accession', inplace=True)

df_refseq = pd.read_csv(REFSEQ_SUMMARY, sep="\t", skiprows=1)
df_refseq.rename(columns={df_refseq.columns[0]: "assembly_accession"}, inplace=True)
df_refseq.set_index('assembly_accession', inplace=True)

df_tree = pd.read_csv(TREE_OF_LIFE)
df_tree.set_index('tax_id', inplace=True)

# dict mapping accession to gff files (in BED format)
with open(f"{BASE_PATH}/slurm/files/gff_bed_files.json") as file:
    gff_files = json.load(file)

def get_gene_data(key):
    # check if gff file exists by looking up in gff_files dictionary
    if key in gff_files:
        gff_bed_file = os.path.join(GFF_BED_PATH, gff_files[key])
        gff = pd.read_csv(gff_bed_file, sep="\t", header=None)
        gff.columns = ['chr', 'start', 'end', 'type']
        # filter by gene and get sum of all gene lengths
        gff = gff[gff['type'] == 'gene']
        gff['length'] = gff['end'] - gff['start']
        gene_content = gff['length'].sum()
        print(f"Gene Content: {gene_content}")
    else:
        gene_content = 0
        print(f"GFF file for {key} does not exist")
    return gene_content
    
def get_species_g4_data(filename, species_data, lock):
    print(f"Getting species data for {filename}")
    try:
        # get species details from assembly summary file
        accession = '_'.join(filename.split('_')[:2])
        if accession.startswith('GCA'):
            df = df_genbank
        elif accession.startswith('GCF'):
            df = df_refseq
        else:
            print(f"Unknown Accession {accession}")
        data = df.loc[accession]
        
        # create new row
        row = {
            'Accession': accession, 
            'Species': data['organism_name'], 
            'Taxa': data['group'], 
            'Genome Size': data['genome_size'], 
            'GC Percentage': data['gc_percent'],
            'Accession': accession.split('_')[1],
            'Gene Count': data['total_gene_count'],
            'Protein Coding Gene Count': data['protein_coding_gene_count'],
            'Non Coding Gene Count': data['non_coding_gene_count'],
            'Gene Content': get_gene_data(accession),
            'Taxid': data['taxid'],
            'Kingdom': df_tree.loc[data['taxid']]['kingdom'],
            'Phylum': df_tree.loc[data['taxid']]['phylum'],
        }
        lock.acquire()
        species_data.append(row)
        lock.release()
        
    except Exception as e:
        print(traceback.format_exc())
        
if __name__ == "__main__":
    # get list of genome files
    with open(FILE_LIST) as f:
        file_list = f.read().splitlines()
    
    # create temp directory if not exists
    if not os.path.exists(f"{BASE_PATH}/temp"):
        os.makedirs(f"{BASE_PATH}/temp")   
    
    # shared variable between processes
    species_data = Manager().list()
    
    # create a lock for the shared variable
    lock = Lock()

    # create processes for each G4 file
    index = 0
    while index < len(file_list):
        processes = []
        for _ in range(NUM_PROC):
            if index == len(file_list):
                break
            process = Process(target=get_species_g4_data, args=(file_list[index], species_data, lock))
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
    df_result = pd.DataFrame(list(species_data))
    df_result.to_csv(RESULT_PATH, index=False)
