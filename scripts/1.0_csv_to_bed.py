import os
import argparse
import pandas as pd
from fuc import pybed

# Directory for temporary bed files
TEMP_PATH = f"temp"

def csv_to_bed(csv_path, file, result_path):
    csv_file = os.path.join(csv_path, file)
    temp_bed_file = os.path.join(TEMP_PATH, file.replace('.csv', '.temp.bed'))
    merged_file = os.path.join(result_path, file.replace('.csv', '.bed'))
    try:
        df = pd.read_csv(csv_file)
        # drop unwanted columns if they exist
        if 'filename' in df.columns:
            df = df.drop('filename', axis=1)
        if 'score' in df.columns: 
            df = df.drop('score', axis=1)
        if 'nbr' in df.columns:
            df = df.drop('nbr', axis=1)
        
        # rename columns in the format required by pybed
        df = df.rename(columns={'chromosome': 'Chromosome', 'start': 'Start', 'end': 'End', 'sequence': 'Sequence', 'length': 'Length'})
        
        print(f"Converting {csv_file}")
        bf = pybed.BedFrame.from_frame(meta=[], data=df)
        bf.to_file(temp_bed_file)

        # sort and merge bed file
        cmd = f"sort -k1,1 -k2,2n {temp_bed_file} | bedtools merge > {merged_file}"
        os.system(cmd)
        
        # remove temporary bed file
        os.remove(temp_bed_file)
        
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert csv to bed files")
    parser.add_argument("csv_path", type=str, help="Path to directory containing CSV files")
    parser.add_argument("file_name", type=str, help="Assembly file name")
    parser.add_argument("bed_path", type=str, help="Path to save bed files")
    cli_args = parser.parse_args()
    
    # directory for resulting merged bed files
    cmd = f"mkdir -p {cli_args.bed_path}"
    os.system(cmd)
    
    csv_to_bed(cli_args.csv_path, cli_args.file_name, cli_args.bed_path)
    