#!/bin/bash
#SBATCH --chdir /storage/group/izg5139/default/akshatha/gquad/data
#SBATCH -o /storage/group/izg5139/default/akshatha/gquad/slurm/logs/extract_out.log
#SBATCH -e /storage/group/izg5139/default/akshatha/gquad/slurm/logs/extract_err.log

source /storage/home/abn5461/miniconda3/bin/activate /storage/home/abn5461/miniconda3/envs/gquad
tar -xf /storage/group/izg5139/default/external/genomes/quadrupia.tar -C /storage/group/izg5139/default/akshatha/gquad/data
ls /storage/group/izg5139/default/akshatha/gquad/data/quadrupia/g4/tar/*.gz |xargs -n1 tar -xzf