#!/bin/bash

j=${1:-1}
latency=${2:-30}

if [[ ! -n "$SSH_CONNECTION" ]];
then
  echo "Local environment detected."
  echo "Initializing bioinformatics genomic compartment bootstrap analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	snakemake --snakefile bedsnake_bootstrap.smk \
            --configfile config_coverage/config.yaml \
            --rerun-incomplete \
	          --rerun-triggers mtime \
            --reason \
            --use-conda \
            --scheduler greedy \
            --keep-going \
            --latency-wait ${latency} \
            --cores $j
else
  echo "SSH Connection detected."
  echo "Initializing bioinformatics genomic compartment bootstrap analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
	snakemake --snakefile bedsnake_bootstrap.smk \
            --configfile config_coverage/config.server.yaml \
            --rerun-incomplete \
	          --rerun-triggers mtime \
            --reason \
            --keep-going \
            --jobs $j \
            --latency-wait ${latency} \
            --cluster-config config_coverage/cluster_bootstrap.yaml \
            --cluster "sbatch -p {cluster.partition} \
                -t {cluster.time} \
                --mem={cluster.mem} \
                --nodes={cluster.nodes} \
                -J {cluster.jobName} \
                -o MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.out \
                -e MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.err"
fi
