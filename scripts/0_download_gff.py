import os
import argparse
import urllib.request

# Download GFF files to this directory
GFF_DIR_PATH = "/storage/group/izg5139/default/external/gff"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert csv to bed files")
    parser.add_argument("path", type=str, help="FTP path for genome assembly")
    cli_args = parser.parse_args()
    
    path = cli_args.path
    path = path.replace("ftp:", "https:")
    try:
        if path.startswith("https:"):
            gff_file = path + "/" + path.split("/")[-1] + "_genomic.gff.gz"
            print(f"Downloading {gff_file}")
            urllib.request.urlretrieve(gff_file, os.path.join(GFF_DIR_PATH, os.path.basename(gff_file)))
            os.system(f"gunzip {os.path.join(GFF_DIR_PATH, os.path.basename(gff_file))}")
            os.remove(os.path.join(GFF_DIR_PATH, os.path.basename(gff_file)))
    except Exception as e:
        print(e)