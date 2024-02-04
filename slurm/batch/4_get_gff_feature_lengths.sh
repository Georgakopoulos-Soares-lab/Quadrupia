#!/bin/bash
#SBATCH --ntasks-per-node=80
#SBATCH -o /storage/group/izg5139/default/akshatha/gquad/slurm/logs/gff_feature_out.log
#SBATCH -e /storage/group/izg5139/default/akshatha/gquad/slurm/logs/gff_feature_err.log

source /storage/home/abn5461/miniconda3/bin/activate /storage/home/abn5461/miniconda3/envs/gquad
python /storage/group/izg5139/default/akshatha/gquad/scripts/4_get_gff_feature_lengths.py