#!/bin/bash
#SBATCH --chdir /storage/group/izg5139/default/akshatha/gquad/data
#SBATCH -o /storage/group/izg5139/default/akshatha/gquad/slurm/logs/jobs/out_%a.log
#SBATCH -e /storage/group/izg5139/default/akshatha/gquad/slurm/logs/jobs/err_%a.log
#SBATCH --array 1-100

# args_file contains the arguments for the script, each line contains the arguments for one job
args_file=/storage/group/izg5139/default/akshatha/gquad/slurm/files/regex_list.txt

# script_file contains the script to run
script_file=/storage/group/izg5139/default/akshatha/gquad/scripts/3.0_get_g4_distribution.py

# path to logs
log_path=/storage/group/izg5139/default/akshatha/gquad/slurm/logs/g4_dist
mkdir -p $log_path

# activate the conda environment
source /storage/home/abn5461/miniconda3/bin/activate /storage/home/abn5461/miniconda3/envs/gquad

# Specify the total number of files using the wc command
total_files=$(cat $args_file | wc -l)
# Specify the number of files to process per job 
# This is equivalent to the number of lines in the file divided by the number of jobs, rounded up
files_per_job=$(echo "($total_files + $SLURM_ARRAY_TASK_COUNT - 1) / $SLURM_ARRAY_TASK_COUNT" | bc)

# Calculate the start and end index for the current array job
start_idx=$(((SLURM_ARRAY_TASK_ID - 1) * files_per_job + 1))
end_idx=$((start_idx + files_per_job - 1))

# Ensure the end index does not exceed the total number of files
if [ $end_idx -gt $total_files ]; then
    end_idx=$total_files
fi

# Process files within the specified range
for i in $(seq $start_idx $end_idx); do
    # get the ith line from the file
    LINE=$(sed -n "$i"p $args_file)
    echo $LINE  
    python $script_file $LINE > $log_path/out_$i.log
done
