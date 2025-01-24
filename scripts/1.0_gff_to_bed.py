import os
import argparse
import traceback
import pandas as pd

# Path to directory containing GFF files
GFF_DIR_PATH = "/storage/group/izg5139/default/external/gff"
# Result path
RESULT_PATH = "gff_bed"

def get_gff_bed_merged(gff_file):
    try:
        file = os.path.basename(gff_file)
        accession = file.replace('.gff', '')
        print(accession, gff_file)
        with open(gff_file) as f:
            lines = f.readlines()
            
        # parse gff file
        # skip metadata lines starting with '#'
        lines = [line.strip() for line in lines if not line.startswith('#')]
        lines = [line.split('\t')[:5] for line in lines]
        df = pd.DataFrame(lines)
        df.columns = ['chr', 'source', 'type', 'start', 'end']
        df = df[['chr', 'start', 'end', 'type']]
        # filter out "region" type
        df = df[df['type'] != 'region']
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df = df.reset_index(drop=True)
        
        # make directory for temp files
        temp_dir = f"temp/{accession}"
        cmd = f"mkdir -p {temp_dir}"
        os.system(cmd)
        
        # GFF files have overlapping regions of the same type
        # merge overlapping regions of the same type separately
        features = df['type'].unique()
        df_dict = {feature: df[df['type'] == feature] for feature in features}
        for feature in features:
            gff_type_bed_file = f"{temp_dir}/{feature}.bed"
            merged_type_bed_file = f"{temp_dir}/{feature}_merged.bed"
            df_dict[feature].to_csv(gff_type_bed_file, sep='\t', index=False, header=False)
            cmd = f"sort -k1,1 -k2,2n {gff_type_bed_file} | bedtools merge > {merged_type_bed_file}"
            os.system(cmd)
            # merged file only retains chr, start, end columns
            # add feature column to merged file
            cmd = f"sed -i 's/$/\t{feature}/' {merged_type_bed_file}"
            os.system(cmd)
            
        # combine all merged bed files into a single sorted bed file
        gff_bed_file = f"{RESULT_PATH}/{accession}.bed"
        cmd = f"cat {temp_dir}/*_merged.bed > {gff_bed_file}"
        os.system(cmd)
        cmd = f"sort -k1,1 -k2,2n {gff_bed_file} -o {gff_bed_file}"
        os.system(cmd)
        
        # GFF is 1-based, convert to 0-based to make it compatible with bed format
        # by subtracting 1 from start index
        # no need to subtract 1 from end because GFF indexes are inclusive
        # while bed indexes are exclusive
        df = pd.read_csv(gff_bed_file, sep='\t', header=None)
        df.columns = ['chr', 'start', 'end', 'feature']
        df['start'] = df['start'].astype(int) - 1
        df.to_csv(gff_bed_file, sep='\t', index=False, header=False)
        
        # remove temp files
        cmd = f"rm -rf {temp_dir}"
        os.system(cmd)
            
    except Exception as e:
        print(traceback.format_exc())
        print(e)
        
if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="Get G4 distribution for each genome")
    parser.add_argument("gff_file", type=str, help="GFF file name")
    cli_args = parser.parse_args()
    
    # directory for resulting merged bed files
    cmd = f"mkdir -p {RESULT_PATH}"
    os.system(cmd)
    
    get_gff_bed_merged(f"{GFF_DIR_PATH}/{cli_args.gff_file}")
