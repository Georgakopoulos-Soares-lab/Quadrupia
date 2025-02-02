# Enrichment TSS/TES Pipeline

__author__ = "Nikol Chantzi"
__email__ = "nmc6088@psu.edu"
__version__ = "1.0.1"

from pathlib import Path
import subprocess
from termcolor import colored
import threading
import csv
import os
import tempfile
import json
import pybedtools
from pybedtools import BedTool
from utils import ProgressTracker
import numpy as np
import pandas as pd
from collections import defaultdict
from mindi.scheduling import MiniBucketScheduler
from mindi.coverage.density import extract_density
from mindi.coverage.pwm_density import strand_evaluators
from mindi.coverage.utils import INTERSECT_FIELDS

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
TOTAL_BUCKETS = int(config['buckets'])
mode = config['mode']
alpha = round(float(config['alpha']), 2)
DESIGN = config['DESIGN']

tempdir = Path(config['tempdir']).resolve()
tempdir.mkdir(exist_ok=True)
pybedtools.helpers.set_tempdir(config['tempdir'])
pybedtools.set_tempdir(config['tempdir'])

def load_bucket(bucket) -> list[str]:
    global TOTAL_BUCKETS
    global out
    with open(f"{out}/schedule_enrichment_{TOTAL_BUCKETS}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

extract_id = lambda assembly: '_'.join(Path(assembly).name.split('_')[:2])

def get_file_ids(design: str) -> dict[str, str]:
    file_ids = {}
    with open(design, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            file_ids.update({row['gff']: (row['accession_id'], row['extraction'])})
    return file_ids

rule all:
    input:
        '%s/%s/enrichment/enrichment_compartments.TSS.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/enrichment_compartments.TES.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/queries_compartments.TSS.%s.csv' % (out, mode, mode),
        '%s/%s/enrichment/queries_compartments.TES.%s.csv' % (out, mode, mode)

rule schedule:
    input:
      DESIGN
    output:
        '%s/schedule_enrichment_%s.json' % (out, TOTAL_BUCKETS)
    run:
        assemblies = []
        with open(input[0], mode='r', encoding='UTF-8') as f:
            reader = csv.DictReader(f, delimiter=",")
            for row in reader:
              assemblies.append(row['gff'])
        print(colored(f"Total assemblies detected: {len(assemblies)}.", "green"))
        if len(assemblies) == 0:
            raise ValueError(f'No assemblies were detected from the design file `{input[0]}`.')
        scheduler = MiniBucketScheduler()
        scheduled_files = scheduler.schedule(assemblies, total_buckets=TOTAL_BUCKETS)
        scheduler.saveas(scheduled_files, output[0])

rule extractEnrichment:
    input:
         DESIGN, 
         '%s/schedule_enrichment_%s.json' % (out, TOTAL_BUCKETS)
    output:
        '%s/%s/enrichment/enrichment_compartments_bucket_{bucket}.%s.enrichment' % (out, mode, mode),
        '%s/%s/enrichment/queries_compartments_bucket_{bucket}.%s.queries' % (out, mode, mode)
    params:
        out=Path(config['out']).resolve(),
        tempdir=Path(config['tempdir']).resolve(),
        compartment=config['compartment'],
        split_category=config['split_category'],
        split_collection=config['split_collection'],
        window_size=int(config["window_size"]),
        sleeping_time=config['log_sleep'],
        mode=config["mode"],
        biotypes=config['biotypes'],
        strand_evaluator=config['strand_evaluator'],
        # bedtools_path=config['bedtools_path'],
    run:
        print(f"Initializing bucket {wildcards.bucket} enrichment extraction process for compartment {params.compartment} and window length {params.window_size}.")
        # >> pybedtools settings
        # pybedtools.set_bedtools_path(path=params.bedtools_path)
        # pybedtools.helpers.set_tempdir("/scratch/nmc6088")
        # pybedtools.helpers.set_bedtools_path(path=params.bedtools_path)
        # << pybedtools settings
      
        # read accessions and map to each gff the accession ID
        accessions = load_bucket(wildcards.bucket)
        file_ids = get_file_ids(input[0])
        # <<
        
        print(colored(f"Splitting coverage process on column {params.split_category}.", "blue"))
        split_category_collection = ["all"]
        if len(params.split_collection) > 0:
            split_category_collection += [int(col) for col in params.split_collection]
        if not params.split_category:
            params.split_category = 'partition'

        total_accessions = len(accessions)

        # >> LOGGING INITIALIZATION
        progress_tracking_log = params.out.joinpath("biologs", f"biolog_tracker_{mode}_{wildcards.bucket}.enrichment.log")
        tracker = ProgressTracker(
                                 total_accessions=total_accessions,
                                 filename=progress_tracking_log,
                                 bucket_id=wildcards.bucket,
                                 sleeping_time=params.sleeping_time)
        logging_thread = threading.Thread(target=tracker.track_progress, daemon=True)
        logging_thread.start()
        # << LOGGING INITIALIZATION

        # >> extraction initializes
        enrichment_table = []
        queries_table = []
        invalid_entries = []
        if len(params.biotypes) == 0:
          biotypes = [None]
        else:
          biotypes = params.biotypes

        for gff_file in accessions:
            tracker.counter += 1
            print(colored(f"Processing accession '{gff_file}'.", "green"))
            accession_id, extraction_file = file_ids[gff_file]

            for biotype in biotypes:
                for attribute in split_category_collection:
                    print(f"Bucket {wildcards.bucket}; Processing {biotype=};{params.split_category}={attribute}...")
                    query_result = extract_density(
                                                extraction=extraction_file,
                                                gff_file=gff_file,
                                                window_size=params.window_size,
                                                compartment=params.compartment,
                                                biotype=biotype,
                                                determine_strand=strand_evaluators[params.strand_evaluator],
                                                attribute_col=params.split_category,
						                        tempdir=tempdir,
                                                attribute=attribute,
                                                enrichment=False,
                                                mode=params.mode
                                                )
                    if query_result is None:
                        print(colored(f"Query failed for {biotype=} and {attribute=} when querying `{gff_file}` GFF.\nThe calculation won't be completed.\nReasons: Empty extraction or empty GFF file.", "red"))
                        invalid_entries.append(extraction_file)
                        continue
                    vector_counts, queries = query_result 
                    # densities table
                    vector_counts.loc[:, params.split_category] = attribute
                    vector_counts.loc[:, "biotype"] = biotype if biotype else '.'
                    vector_counts.loc[:, "#assembly_accession"] = accession_id
                    # queries table
                    queries.index.name = "site"
                    queries = queries.reset_index()
                    queries.loc[:, "#assembly_accession"] = accession_id
                    queries.loc[:, params.split_category] = attribute
                    queries.loc[:, "biotype"] = biotype if biotype else '.'
                    # column rearrangement
                    prev_cols = set(vector_counts.columns)
                    vector_counts = vector_counts[["#assembly_accession",
                                                   "site",
                                                   "biotype",
                                                   params.split_category,
                                                   ] + list(range(-params.window_size, params.window_size+1))]
                    if params.mode == "template":
                        vector_counts = vector_counts.reset_index()\
                                                     .rename(columns={"index": "template|non_template"})
                    enrichment_table.append(vector_counts)
                    print(vector_counts)
                    queries_table.append(queries)
        # << extraction finished
        
        # >> save results
        if len(queries_table) > 0:
            queries_table = pd.concat(queries_table, axis=0)
        else:
            queries_table = pd.DataFrame([], columns=["#assembly_accession",
                                                     "site",
                                                     "total_targets",
                                                     "matched_targets",
                                                     "matched_targets_non_template",
                                                     "matched_targets_template",
                                                     "target_perc",
                                                     "total_queries",
                                                     "matched_queries",
                                                     "query_perc",
                                                     "partition",
                                                     "biotype"])
        if len(enrichment_table) > 0:
            enrichment_table = pd.concat(enrichment_table, axis=0)
        else:
            print(colored(f"No files detected for bucket {wildcards.bucket}.", "red"))
            enrichment_table = pd.DataFrame([], columns=["#assembly_accession",
                                                         "site",
                                                         "biotype",
                                                         params.split_category,
                                                         "nucleotide"] + list(range(-params.window_size, params.window_size+1))
                                            )
        pybedtools.helpers.cleanup(remove_all=False)
        # save enrichment table
        enrichment_table.set_index("#assembly_accession", inplace=True)
        enrichment_table.to_csv(output[0], sep=",", index=True, mode="w")

        # save queries table 
        queries_table.set_index("#assembly_accession", inplace=True)
        queries_table.to_csv(output[1], sep=",", index=True, mode="w")
        # << save results finished

rule reduceEnrichment:
    input:
        expand([
               '%s/%s/enrichment/enrichment_compartments_bucket_{bucket}.%s.enrichment' % (out, mode, mode),
               '%s/%s/enrichment/queries_compartments_bucket_{bucket}.%s.queries' % (out, mode, mode)
               ],
               bucket=range(TOTAL_BUCKETS))
    output:
        '%s/%s/enrichment/enrichment_compartments.TSS.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/enrichment_compartments.TES.%s.parquet' % (out, mode, mode),
        '%s/%s/enrichment/queries_compartments.TSS.%s.csv' % (out, mode, mode),
        '%s/%s/enrichment/queries_compartments.TES.%s.csv' % (out, mode, mode),
    run:
        # enrichment
        enrichment_table = []
        for bucket in range(TOTAL_BUCKETS):
            coverage_file = f"{out}/{mode}/enrichment/enrichment_compartments_bucket_{bucket}.{mode}.enrichment"
            coverage_df = pd.read_csv(coverage_file)
            enrichment_table.append(coverage_df)
        enrichment_table = pd.concat(enrichment_table, axis=0)
        enrichment_table.set_index("#assembly_accession", inplace=True)
        enrichment_table_TSS = enrichment_table[enrichment_table['site'] == 'tss'].drop(columns=['site'])
        enrichment_table_TES = enrichment_table[enrichment_table['site'] == 'tes'].drop(columns=['site'])
        # save enrichment table
        enrichment_table_TSS.to_parquet(output[0], engine="fastparquet")
        enrichment_table_TES.to_parquet(output[1], engine="fastparquet")
        # queries
        queries_table = []
        for bucket in range(TOTAL_BUCKETS):
            coverage_file = f"{out}/{mode}/enrichment/queries_compartments_bucket_{bucket}.{mode}.queries"
            coverage_df = pd.read_csv(coverage_file)
            queries_table.append(coverage_df)
        queries_table = pd.concat(queries_table, axis=0)
        queries_table.set_index("#assembly_accession", inplace=True)
        queries_table_TSS = queries_table[queries_table['site'] == 'tss'].drop(columns=['site'])
        queries_table_TES = queries_table[queries_table['site'] == 'tes'].drop(columns=['site'])
        # save query table
        queries_table_TSS.to_csv(output[2], sep=",", index=True, mode="w")
        queries_table_TES.to_csv(output[3], sep=",", index=True, mode="w")
