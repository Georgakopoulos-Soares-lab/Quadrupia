import os

# Base path
BASE_PATH = "/storage/group/izg5139/default/akshatha/gquad"
# list of files
FILE_LIST = "/storage/group/izg5139/default/akshatha/gquad/slurm/files/g4_list.txt"

def get_overlapping_seq(filename):
    g4_bed_file = f"{BASE_PATH}/data/g4hunter_bed/{filename}"
    reg_bed_file = f"{BASE_PATH}/data/regex_bed/{filename}"
    int_bed_file = f"{BASE_PATH}/temp/{filename}"

    g4_count = 0
    reg_count = 0
    int_count = 0
    
    try:
        # count the number of lines in each file for genomic files
        cmd = f"wc -l {g4_bed_file}"
        g4_count = int(os.popen(cmd).read().split()[0])
        
        cmd = f"wc -l {reg_bed_file}"
        reg_count = int(os.popen(cmd).read().split()[0])
        
        # find overlapping sequences between g4hunter and regex
        cmd = f"bedtools intersect -a {g4_bed_file} -b {reg_bed_file} -u -f 0.50 -r > {int_bed_file}"
        os.system(cmd)
        cmd = f"wc -l {int_bed_file}"
        int_count = int(os.popen(cmd).read().split()[0])
    
        os.remove(int_bed_file)
        
    except Exception as e:
        print(e)
        
    print(g4_count, reg_count, int_count)
    counts = {
        "g4hunter": g4_count,
        "regex": reg_count,
        "intersect": int_count
    }
    
    return counts

if __name__ == '__main__':
    c = {
        "g4hunter": 0,
        "regex": 0,
        "intersect": 0
    }
    
    with open(FILE_LIST) as f:
        file_list = f.readlines()
        file_list = [file.strip() for file in file_list]
        
    for file in file_list:
        counts = get_overlapping_seq(file.replace('csv', 'bed'))
        c['g4hunter'] += counts['g4hunter']
        c['regex'] += counts['regex']
        c['intersect'] += counts['intersect']
        
    with open(f"{BASE_PATH}/results/counts.txt", 'w') as f:
        f.write(f"G4Hunter: {c['g4hunter']}\n")
        f.write(f"Regex: {c['regex']}\n")
        f.write(f"Intersect: {c['intersect']}")
