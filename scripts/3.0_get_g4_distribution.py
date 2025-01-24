import os
import json
import argparse
import traceback
import pandas as pd

# G4 Files
G4_PATH = "data/g4hunter_bed"
REGEX_PATH = "data/regex_bed"
# GFF files(in BED format)
GFF_BED_PATH = "data/gff_bed"

def get_g4_data(g4_file):
    file = os.path.basename(g4_file)
    g4_data = {
        "accession": '_'.join(file.split('_')[:2]),
        "g4_count": 0,
        "total_g4_length": 0
    }
    try:
        # get no. of G4s in the file
        g4_data['g4_count'] = int(os.popen(f"wc -l {g4_file}").read().split()[0])
        
        # get total length of G4s
        with open(g4_file) as f:
            data = f.readlines()
            data = [line.split('\t') for line in data]
        data = [[int(line[1]), int(line[2])] for line in data]
        g4_data['total_g4_length'] = sum([end - start for start, end in data])
    except Exception as e:
        print("Error: ", e)
        print(traceback.format_exc())
    finally:
        return g4_data
    

def get_g4_distribution(g4_file, g4_data):
    print(f"Getting G4 distribution for file {g4_file}")
    file = os.path.basename(g4_file)
    
    try:
        # get overlap between G4 and GFF using bedtools intersect
        # keep all G4 data
        intersect_file = f"{BASE_PATH}/temp/{file}_intersect.bed"
        cmd = f"bedtools intersect -a {g4_file} -b {gff_bed_file} -wao > {intersect_file}"
        os.system(cmd)
        
        # clean up intersect file
        with open(intersect_file) as f:
            lines = f.readlines()
            lines = [line.split('\t')[:8] for line in lines]
        with open(intersect_file, 'w') as f:
            for line in lines:
                f.write('\t'.join(line) + '\n')
        
        # read intersected data into a dataframe
        df = pd.read_csv(intersect_file, sep='\t', header=None)
        df.columns = ['g4_chr', 'g4_start', 'g4_end', 'gff_chr', 'gff_start', 'gff_end', 'feature', 'overlap_length']

        # NaN in feature column indicate G4s that do not overlap with any GFF feature
        # hence these lie in the intergenic region
        df['feature'] = df['feature'].fillna('intergenic')
        
        # convert start and end to int
        df['g4_start'] = df['g4_start'].astype(int)
        df['g4_end'] = df['g4_end'].astype(int)
        df['gff_start'] = df['gff_start'].astype(int)
        df['gff_end'] = df['gff_end'].astype(int)
        
        # overlap length of '.' indicates no overlap, i.e. this entire G4 falls in intergenic region
        # hence overlap length is the length of the G4
        df['overlap_length'] = df.apply(lambda x: x['g4_end'] - x['g4_start'] if x['overlap_length'] == '.' else x['overlap_length'], axis=1)
        df['overlap_length'] = df['overlap_length'].astype(int)
        
        # get sum of overlap lengths for each feature
        df = df.groupby(['feature']).agg({'overlap_length': 'sum'}).reset_index()
        total_g4_dist = df.set_index('feature').T.to_dict('records')[0]
        
        # add accession, total G4 count and length to the dictionary
        total_g4_dist['accession'] = g4_data['accession']
        total_g4_dist['g4_count'] = g4_data['g4_count']
        total_g4_dist['total_g4_length'] = g4_data['total_g4_length']
        
        # remove intersect file
        os.remove(intersect_file)
        
    except Exception as e:
        print("Error: ", e)
        print(traceback.format_exc())
        total_g4_dist = g4_data.copy()
        
    finally:
        return total_g4_dist
        
if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="Get G4 distribution for each genome")
    parser.add_argument("g4_file", type=str, help="G4 file name")
    cli_args = parser.parse_args()
    
    filename = cli_args.g4_file
    filename = filename.replace('csv', 'bed')
    
    # get no. of G4s and total G4 content from G4 BED files
    g4hunter_data = get_g4_data(f"{G4_PATH}/{filename}")
    regex_data = get_g4_data(f"{REGEX_PATH}/{filename}")
    
    # check if gff file exists by looking up in gff_files dictionary
    accession = '_'.join(filename.split('_')[:2])
    with open(f"{BASE_PATH}/slurm/files/gff_bed_files.json") as gff_map_file:
        gff_files = json.load(gff_map_file)
    if accession in gff_files:
        gff_bed_file = f"{GFF_BED_PATH}/{gff_files[accession]}"
        
        # add G4 distribution data to G4 data dictionary
        total_g4_dist ={
            "g4hunter": get_g4_distribution(f"{G4_PATH}/{filename}", g4hunter_data),
            "regex": get_g4_distribution(f"{REGEX_PATH}/{filename}", regex_data)
        }

    else:
        total_g4_dist = {
            "g4hunter": g4hunter_data,
            "regex": regex_data
        }
        print(f"GFF file not found for {filename}")
        
    # save dictionary to file
    cmd = f"mkdir -p {BASE_PATH}/data/g4_distribution"
    os.system(cmd)
    result_file = f"{BASE_PATH}/data/g4_distribution/{accession}.json"
    with open(result_file, 'w') as file:
        json.dump(total_g4_dist, file)
