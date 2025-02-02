from tqdm import tqdm
from termcolor import colored
import numpy as np
import pandas as pd
import polars as pl
from pathlib import Path

if __name__ == "__main__":

    # Script that calculates densities of motifs from extracted genomes for various taxonomic ranks
    import argparse

    parser = argparse.ArgumentParser(description="""Calculates densities for domains & phylums.""")
    parser.add_argument("--dest", "-d", type=str, default="density")
    parser.add_argument("--source", "-s", type=str, required=True)
    parser.add_argument("--empty", "-e", type=str, required=True)
    parser.add_argument("--seq_col", "-sc", type=str, default="sequenceOfArm")
    parser.add_argument("--step_threshold", "-st", type=float, default=0.05)
    parser.add_argument("--merged", type=int, default=0)
    parser.add_argument("--calculate_gc", type=int, default=1)
    parser.add_argument("--mode", "-m", type=str,
                        required=True,
                        choices=["STR", "IR", "MR", "DR", "Z", "G4", "H-DNA"])

    args = parser.parse_args()
    source = args.source
    empty = args.empty
    mode = args.mode
    merged = args.merged
    calculate_gc = args.calculate_gc
    seq_col = args.seq_col
    step_threshold = args.step_threshold
    assert step_threshold >= 0 and step_threshold < 1.0, f"GC Content step threshold must lie within [0, 1) range, but `{step_threshold}` was provided."
    dest = Path(args.dest).resolve().joinpath("density")
    dest.mkdir(exist_ok=True, parents=True)

    # load empty assemblies
    columns = ["#assembly_accession",
               "superkingdom", 
               "kingdom", 
               "phylum", 
               "order", 
               "class", 
               "family", 
               "genome_size", 
               "group",
               ]
    df_empty = pl.read_csv(empty, separator="\t")\
                    .filter(pl.col(mode) == 1)\
                    .select(columns)
    empty_set = set(df_empty["#assembly_accession"])
    # load not-empty assemblies
    df = pl.read_parquet(source)
    if "sequenceLength" not in df.columns:
        df = df.with_columns(
                (pl.col("end") - pl.col("start")).alias("sequenceLength")
                )

    def calculate_density(df: pl.DataFrame, df_empty) -> None:
        empty_set = set(df_empty["#assembly_accession"])
        df_grouped = df.group_by("#assembly_accession", maintain_order=True)\
            .agg([
                    pl.col("superkingdom").first(),
                    pl.col("kingdom").first(),
                    pl.col("phylum").first(),
                    pl.col("order").first(),
                    pl.col("class").first(),
                    pl.col("family").first(),
                    pl.col("genome_size").first(),
                    # pl.col("gc_percent").first(),
                    pl.col("group").first(),
                    pl.col("start").count().alias("totalCounts"),
                    pl.col("sequenceLength").sum().alias("total_bp")
                    ])\
            .with_columns(
                    (pl.col("totalCounts") * 1e6 / pl.col("genome_size"))\
                            .round(2)\
                            .alias("density")
                    )\
            .sort("density", descending=True)
        non_empty_set = set(df_grouped["#assembly_accession"])
        assert len(non_empty_set.intersection(empty_set)) == 0, "Empty assemblies cannot overlap with non-empty assemblies."

        # define constants: 
        # - total number of assemblies 
        all_assemblies = set(df_grouped["#assembly_accession"].unique())
        total_assemblies = len(all_assemblies)
    
        # this table contains all assemblies!
        prev = df_grouped.shape[0]
        df_grouped = pl.concat([df_grouped, df_empty], how="diagonal")\
                    .with_columns(
                            pl.col("density").fill_null(0.0)
                            )
        assert df_grouped.shape[0] == prev + df_empty.shape[0], "Failed to merge properly."

        # save table 
        df_grouped.write_csv(f"{dest}/density_data_merged_{merged}_{mode}.txt", separator="\t")
        ranks = ["superkingdom",
                 "kingdom",
                 "phylum",
                 "order",
                 "class",
                 "family"]

        print(colored("Calculating densities for taxonomic ranks.", "green"))
        for rank in ranks:
            df_grouped_rank = df_grouped\
                        .filter(
                                ~pl.all_horizontal(pl.col(rank).is_null())
                        )\
                        .group_by(rank, maintain_order=True)\
                    .agg([
                            pl.col("density").mean().round(2).alias("avg_density"),
                            pl.col("density").std().round(2).alias("std_density"),
                            pl.col("density").median().round(2).alias("median_density"),
                            pl.col("totalCounts").sum(),
                            ])\
                    .sort("avg_density", descending=True)
            df_grouped_rank.write_csv(f"{dest}/density_data_merged_{merged}_{mode}_{rank}.txt", separator="\t")
        print(colored("Density calculation for taxonomic ranks has been succesfully completed.", "green"))
        return

    df_grouped = df.group_by("#assembly_accession", maintain_order=True)\
        .agg([
                pl.col("superkingdom").first(),
                pl.col("kingdom").first(),
                pl.col("phylum").first(),
                pl.col("order").first(),
                pl.col("class").first(),
                pl.col("family").first(),
                pl.col("genome_size").first(),
                # pl.col("gc_percent").first(),
                pl.col("group").first(),
                pl.col("start").count().alias("totalCounts"),
                pl.col("sequenceLength").sum().alias("total_bp")
                ])\
        .with_columns(
                (pl.col("totalCounts") * 1e6 / pl.col("genome_size"))\
                        .round(2)\
                        .alias("density")
                )\
        .sort("density", descending=True)
    non_empty_set = set(df_grouped["#assembly_accession"])
    assert len(non_empty_set.intersection(empty_set)) == 0, "Empty assemblies cannot overlap with non-empty assemblies."
    
    # define constants: 
    # - total number of assemblies 
    all_assemblies = set(df_grouped["#assembly_accession"].unique())
    total_assemblies = len(all_assemblies)
    
    # this table contains all assemblies!
    prev = df_grouped.shape[0]
    df_grouped = pl.concat([df_grouped, df_empty], how="diagonal")\
                    .with_columns(
                            pl.col("density").fill_null(0.0)
                            )
    assert df_grouped.shape[0] == prev + df_empty.shape[0], "Failed to merge properly."

    # save table 
    df_grouped.write_csv(f"{dest}/density_data_merged_{merged}_{mode}.txt", separator="\t")
    ranks = ["superkingdom",
             "kingdom",
             "phylum",
             "order",
             "class",
             "family"]

    print(colored("Calculating densities for taxonomic ranks.", "green"))
    for rank in ranks:
        df_grouped_rank = df_grouped\
                        .filter(
                                ~pl.all_horizontal(pl.col(rank).is_null())
                        )\
                        .group_by(rank, maintain_order=True)\
                    .agg([
                            pl.col("density").mean().round(2).alias("avg_density"),
                            pl.col("density").std().round(2).alias("std_density"),
                            pl.col("density").median().round(2).alias("median_density"),
                            pl.col("totalCounts").sum(),
                            ])\
                    .sort("avg_density", descending=True)
        df_grouped_rank.write_csv(f"{dest}/density_data_merged_{merged}_{mode}_{rank}.txt", separator="\t")

    # calculate rolling density bound on GC threshold
    # calculate GC content
    if calculate_gc == 0:
        import sys
        sys.exit(0)

    df = df.with_columns(
                    pl.col(seq_col).map_elements(lambda seq: (seq.count("g") + seq.count("c"))/len(seq),
                                                 return_dtype=float
                                                 )
                                   .alias("gc_content")
                    )
    print(colored("Density calculation for taxonomic ranks has been succesfully completed.", "green"))

    # load empty and non-empty assemblies for all arms
    empty_total = {}
    thresholds = np.arange(0.0, 1.0, step_threshold)
    df_thresh = df
    df_grouped_rank_all = []
    print(colored("Calculating densities for taxonomic ranks for various GC content thresholds on sequence motif.", "green"))
    for threshold in tqdm(thresholds, leave=True, position=0):
        df_thresh = df_thresh.filter(pl.col("gc_content") >= threshold)
        df_grouped_remaining = df_thresh.group_by("#assembly_accession", 
                                                  maintain_order=True)\
                                    .agg(
                                            pl.col("superkingdom").first(),
                                            pl.col("kingdom").first(),
                                            pl.col("phylum").first(),
                                            pl.col("genome_size").first(),
                                            # pl.col("gc_percent").first(),
                                            pl.col("group").first(),
                                            pl.col("start").count().alias("totalCounts"),
                                            pl.col("sequenceLength").sum().alias("total_bp")
                                    )\
                                .with_columns(
                                        (pl.col("totalCounts") * 1e6 / pl.col("genome_size"))\
                                                .round(2)\
                                                .alias("density")
                                        )\
                                .sort("density", descending=True)

        # determine empty accessions with the new threshold level
        # empty is everything not found in df_grouped
        empty_assemblies = list(all_assemblies - set(df_grouped_remaining["#assembly_accession"]))
        df_empty = pl.DataFrame({"#assembly_accession": empty_assemblies})
        if df_empty.shape[0] > 0:
            df_empty = df_empty.join(df_grouped\
                                .select(["#assembly_accession",
                                         "genome_size",
                                         "superkingdom",
                                         "kingdom",
                                         "phylum",
                                         "class",
                                         "order",
                                         "family"
                                         ]), 
                              on="#assembly_accession")
        # ensure all assemblies have been succesfully mapped to their metadata
        assert df_empty.shape[0] == len(empty_assemblies), f"Expected empty assemblies: {df_empty.shape[0]} vs. total {len(empty_assemblies)}."
        empty_total[threshold] = len(empty_assemblies)

        # join empty and non-empty to estimate motif density 
        if df_empty.shape[0] > 0:
            df_grouped_remaining = pl.concat([df_grouped_remaining, df_empty], how="diagonal")\
                                        .with_columns(
                                                pl.col("density").fill_null(0.0)
                                            )
        assert df_grouped_remaining.shape[0] == total_assemblies, f"Expected empty assemblies: {df_grouped_remaining.shape[0]} vs. total {total_assemblies}."
        # estimate density for each rank on GC threshold
        # from superkingdom, kingdom down to phylum level suffices?
        for rank in ["superkingdom", "kingdom"]:
            df_grouped_rank = df_grouped_remaining\
                            .filter(
                                    ~pl.all_horizontal(pl.col(rank).is_null())
                            )\
                            .group_by(rank, maintain_order=True)\
                            .agg([
                                    pl.col("density").mean().round(2).alias("avg_density"),
                                    pl.col("density").std().round(2).alias("std_density"),
                                    pl.col("density").median().round(2).alias("median_density"),
                                    pl.col("totalCounts").sum()
                                    ])\
                            .sort("avg_density", descending=True)\
                            .rename({rank: "taxonomic_rank"})\
                            .with_columns(
                                    pl.lit(threshold).alias("gc_threshold")
                            ).to_pandas()
            df_grouped_rank_all.append(df_grouped_rank)

    df_grouped_rank_all = pd.concat(df_grouped_rank_all)
    df_grouped_rank_all.to_csv(f"{dest}/density_data_merged_{merged}_taxonomic_ranks_varying_gc_threshold_{mode}.txt", 
                                sep="\t",
                                mode="w",
                                index=False,
                                header=True)

    # save empty total per threshold
    with open(f"{dest}/empty_assemblies_varying_gc_threshold_{mode}.txt", mode="w", encoding="utf-8") as f:
        for gc_threshold, total_empty in empty_total.items():
            f.write(f"{gc_threshold}\t{total_empty}\n")
