#!/bin/bash
#SBATCH --ntasks-per-node=80
#SBATCH -o /storage/group/izg5139/default/akshatha/gquad/slurm/logs/species_out.log
#SBATCH -e /storage/group/izg5139/default/akshatha/gquad/slurm/logs/species_err.log

source /storage/home/abn5461/miniconda3/bin/activate /storage/home/abn5461/miniconda3/envs/gquad
python /storage/group/izg5139/default/akshatha/gquad/scripts/2_get_species_data.py