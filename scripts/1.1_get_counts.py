import os

# list of files
FILE_LIST = "/storage/group/izg5139/default/akshatha/gquad/slurm/files/g4_list.txt"

def get_overlapping_seq(filename):
    g4_bed_file = f"data/g4hunter_bed/{filename}"
    reg_bed_file = f"data/regex_bed/{filename}"
    int_bed_file1 = f"temp/1_{filename}"
    int_bed_file2 = f"temp/2_{filename}"

    g4_count = 0
    reg_count = 0
    int_count1 = 0
    int_count2 = 0
    
    try:
        # count the number of lines in each file for genomic files
        cmd = f"wc -l {g4_bed_file}"
        g4_count = int(os.popen(cmd).read().split()[0])
        
        cmd = f"wc -l {reg_bed_file}"
        reg_count = int(os.popen(cmd).read().split()[0])
        
        # find overlapping sequences between g4hunter and regex with atleast 1bp overlap
        # -a : file A
        # -b : file B
        # -u : write the original A entry once if any overlap with B
        cmd = f"bedtools intersect -a {g4_bed_file} -b {reg_bed_file} -u > {int_bed_file1}"
        os.system(cmd)
        cmd = f"wc -l {int_bed_file1}"
        int_count1 = int(os.popen(cmd).read().split()[0])
        os.remove(int_bed_file1)
        
        # find overlapping sequences between g4hunter and regex
        # -a : file A
        # -b : file B
        # -u : write the original A entry once if any overlap with B
        # -f : require that the overlap be at least 50% of the A feature
        # -F : require that the overlap be at least 50% of the B feature
        # -e : require that the minimum fraction be satistied for either A or B
        # that is, either A or B must have at least 50% overlap
        cmd = f"bedtools intersect -a {g4_bed_file} -b {reg_bed_file} -u -f 0.50 -F 0.50 -e > {int_bed_file2}"
        os.system(cmd)
        cmd = f"wc -l {int_bed_file2}"
        int_count2 = int(os.popen(cmd).read().split()[0])
        os.remove(int_bed_file2)
        
    except Exception as e:
        print(e)
        
    print(g4_count, reg_count, int_count1, int_count2)
    counts = {
        "g4hunter": g4_count,
        "regex": reg_count,
        "intersect_min1": int_count1,
        "intersect_min50": int_count2
    }
    
    return counts

if __name__ == '__main__':
    c = {
        "g4hunter": 0,
        "regex": 0,
        "intersect_min1": 0,
        "intersect_min50": 0
    }
    
    with open(FILE_LIST) as f:
        file_list = f.readlines()
        file_list = [file.strip() for file in file_list]
        
    # create temp directory
    if not os.path.exists("temp"):
        os.makedirs("temp")
        
    for file in file_list:
        counts = get_overlapping_seq(file.replace('csv', 'bed'))
        c['g4hunter'] += counts['g4hunter']
        c['regex'] += counts['regex']
        c['intersect_min1'] += counts['intersect_min1']
        c['intersect_min50'] += counts['intersect_min50']
        
    with open("results/counts.txt", 'w') as f:
        f.write(f"G4Hunter: {c['g4hunter']}\n")
        f.write(f"Regex: {c['regex']}\n")
        f.write(f"Intersect (min 1bp): {c['intersect_min1']}\n")
        f.write(f"Intersect (min 50%): {c['intersect_min50']}")
