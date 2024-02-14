#!/bin/bash
#SBATCH --chdir /storage/group/izg5139/default/akshatha/gquad
#SBATCH -o /storage/group/izg5139/default/akshatha/gquad/slurm/logs/count_out.log
#SBATCH -e /storage/group/izg5139/default/akshatha/gquad/slurm/logs/count_err.log

source /storage/home/abn5461/miniconda3/bin/activate /storage/home/abn5461/miniconda3/envs/gquad
python /storage/group/izg5139/default/akshatha/gquad/scripts/1.1_get_counts.py
