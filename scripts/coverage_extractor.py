# GFF Utils + Coverage Workflow

# Yet another GFF parser built on Polars

from pathlib import Path
import time
import argparse
from tqdm import tqdm
import csv
import gzip
import os
from termcolor import colored
import json
import logging
import polars as pl
import pandas as pd
import threading
import pysam
import attr
from attr import field
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Optional, ClassVar, Iterator
from pybedtools import BedTool


@attr.s(slots=True, kw_only=True)
class GFFExtractor:

    compartments: Optional[list[str]] = field(factory=lambda : ["gene", "CDS", "exon", "region", "pseudogene", "genes", "pseudogenes"])
    names_mapping: dict[str, str] = field(factory=lambda : {"gene": "Gene",
                                                            "exon": "Exon",
                                                            "region": "Genome",
                                                            "intron": "Intron",
                                                            "pseudogene": "Pseudogene",
                                                            })
    promoter_kb: int = field(default=1_000, validator=attr.validators.instance_of(int), converter=int)
    terminator_kb: int = field(default=1_000, validator=attr.validators.instance_of(int), converter=int)
    delim: str = field(default="*", converter=str)
    GFF_FIELDS: ClassVar[list[str]] = ["seqID", "source", "compartment", "start", "end", "score", "strand", "phase", "attributes"]
    COVERAGE_FIELDS: ClassVar[list[str]] = ["totalHits", "overlappingBp", "compartmentLength", "coverage"]

    def __attrs_post_init__(self) -> None:
        pass

    @staticmethod
    def parse_attributes(attributes: str) -> dict[str, str]:
        attributes = attributes.split(";")
        attributes_map = {}
        for attr in attributes:
            key, value = attr.split("=")
            attributes_map[key] = value
        return attributes_map

    @staticmethod 
    def parse_biotype(attributes: str) -> str:
        attributes = GFFExtractor.parse_attributes(attributes)
        return attributes.get("biotype", attributes.get("gene_biotype", "."))

    @staticmethod 
    def parse_ID(attributes: str) -> str:
        attributes = GFFExtractor.parse_attributes(attributes)
        return attributes.get("ID", ".")
    
    @staticmethod
    def parse_gene_biotype(gff_table: pl.DataFrame) -> pl.DataFrame:
        return gff_table.with_columns(
                                pl.when(pl.col("compartment") == "gene")
                                .then(
                                    pl.when(pl.col("attributes").str.extract(r"gene_biotype=([^;]+)", 1) == "protein_coding")
                                    .then(pl.lit("protein_coding"))
                                    .otherwise(pl.lit("non_coding"))
                                    )
                                .otherwise(pl.lit("."))
                                .alias("biotype")
                            )

    @staticmethod
    def parse_compartment_ID(gff_table: pl.DataFrame) -> pl.DataFrame:
            return gff_table.with_columns(
                            pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("compartmentID")
                        )

    @staticmethod
    def parse_parent_ID(gff_table: pl.DataFrame) -> pl.DataFrame:
            return gff_table.with_columns(
                            pl.col("attributes").str.extract(r"ParentID=([^;]+)", 1).alias("parentID")
                        )

    def read_gff(self, gff_file: str,
                        change_names: bool = True,
                        parse_biotype: bool = True,
                        parse_ID: bool = False,
                        parse_parentID: bool = False,
                        replace_pseudogene_with_gene: bool = True,
                        parse_promoters: bool = False,
                        parse_terminators: bool = False,
                 ) -> pl.DataFrame:
        gff_table = pl.read_csv(gff_file, 
                                separator="\t",
                                has_header=False,
                                comment_prefix="#",
                                new_columns=GFFExtractor.GFF_FIELDS
                                )\
                        .with_columns(
                                        # 1-base --> 0-base
                                        pl.col("start") - 1,
                                        pl.col("compartment").replace({
                                                                       "genes": "gene", 
                                                                       "pseudogenes": "pseudogene"
                                                                       })
                            )
        gff_table = gff_table.join(
                        gff_table.filter(pl.col("compartment") == "region")\
                                 .select(["seqID", "start", "end"])\
                                 .rename(mapping={
                                            "start": "region_start",
                                            "end": "region_end",
                                            }
                                         ),
                        on="seqID",
                        how="left",
                        )
        if self.compartments:
            compartments = set(self.compartments)
            gff_table = gff_table.filter(pl.col("compartment").is_in(compartments))
        if replace_pseudogene_with_gene:
            gff_table = gff_table.with_columns(
                                pl.col("compartment").replace("pseudogene", "gene").alias("compartment")
                                )
        if parse_biotype:
            gff_table = GFFExtractor.parse_gene_biotype(gff_table)
        if parse_ID:
            gff_table = GFFExtractor.parse_compartment_ID(gff_table)
        if parse_parentID:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
        if change_names:
            gff_table = gff_table.with_columns(
                                compartment=pl.col("compartment").replace(self.names_mapping)
                                    )
        if parse_promoters:
            gff_table = pl.concat([
                            gff_table,
                            self.parse_promoters(gff_table, 
                                                 promoter_kb=self.promoter_kb)
                            ])
        if parse_terminators:
            gff_table = pl.concat([
                            gff_table,
                            self.parse_terminators(gff_table, 
                                                   terminator_kb=self.terminator_kb)
                            ])
        return gff_table

    def parse_terminators(self, gff_table: pl.DataFrame, terminator_kb: int = 1_000, use_names: bool = False) -> pl.DataFrame:
        selection_cols = ["seqID", "terminator_start", "terminator_end", "strand", "region_start", "region_end"]
        if use_names:
            selection_cols.insert(5, "compartmentID")
        return gff_table.filter(pl.col("end") <= pl.col("region_end"), pl.col("strand") == "+" | pl.col("strand") == "-")\
                        .with_columns(
                                    [
                                    pl.when(pl.col("strand") == "+")
                                        .then(pl.col("end"))
                                        .otherwise(pl.max_horizontal(0, pl.col("start") - terminator_kb))
                                        .alias("terminator_start"),
                                    pl.when(pl.col("strand") == "+")
                                        .then(pl.min_horizontal(pl.col("end") + terminator_kb, pl.col("region_end")))
                                        .otherwise(pl.col("start"))
                                        .alias("terminator_end")
                                    ]
                                )\
                        .filter(pl.col("terminator_start") < pl.col("terminator_end"))\
                        .select(selection_cols)\
                        .rename({
                                 "terminator_start": "start",
                                 "terminator_end": "end"
                                 })\
                        .with_columns(
                                pl.lit("Terminator").alias("compartment")
                                )

    def parse_promoters(self, gff_table: pl.DataFrame, promoter_kb: int = 1_000, use_names: bool = False) -> pl.DataFrame:
        selection_cols = ["seqID", "promoter_start", "promoter_end", "strand", "region_start", "region_end"]
        if use_names:
            selection_cols.insert(5, "compartmentID")
        return gff_table.filter(pl.col("end") <= pl.col("region_end"), pl.col("strand") == "+" | pl.col("strand") == "-")\
                             .with_columns(
                                            [
                                            pl.when(pl.col("strand") == "+")
                                                .then(pl.max_horizontal(0, pl.col("start") - promoter_kb))
                                                .otherwise(pl.col("end"))
                                                .alias("promoter_start"),
                                            pl.when(pl.col("strand") == "+")
                                                .then(pl.col("start"))
                                                .otherwise(pl.min_horizontal(pl.col("end") + promoter_kb, pl.col("region_start")))
                                                .alias("promoter_end")
                                                ]
                                            )\
                        .filter(pl.col("promoter_start") < pl.col("promoter_end"))\
                        .select(selection_cols)\
                        .rename({
                                 "promoter_start": "start",
                                 "promoter_end": "end"
                                 })\
                        .with_columns(
                                pl.lit("Promoter").alias("compartment")
                            )

    def merge(self, gff_table: pl.DataFrame, faidx: Optional[str] = None) -> pl.DataFrame:
        unique_compartments = set(gff_table["compartment"])
        merged_gff: list[pl.DataFrame] = []
        for compartment in unique_compartments:
            gff_table_bed = BedTool.from_dataframe(
                                        gff_table.filter(pl.col("compartment") == compartment)
                                        .select(["seqID", "start", "end"])
                                        .to_pandas()
                                        )
            if faidx is not None:
                gff_table_bed = gff_table_bed.sort(faidx=faidx)
            else:
                gff_table_bed = gff_table_bed.sort()
            gff_temp = pl.read_csv(
                             gff_table_bed.merge(c="3", 
                                                 o="count", 
                                                 delim=self.delim).fn,
                             has_header=False,
                             new_columns=["seqID", "start", "end", "counts"],
                             separator="\t",
                            )\
                            .with_columns(
                                pl.lit(compartment).alias("compartment")
                            )
            merged_gff.append(gff_temp)
        merged_gff = pl.concat(merged_gff)
        return merged_gff

    def parse_coverage(self, gff_table: pl.DataFrame, 
                            extraction_table: pl.DataFrame,
                            group: bool = True,
                            faidx: Optional[str] = None) -> pl.DataFrame:
        if isinstance(gff_table, pl.DataFrame):
            gff_table = gff_table.to_pandas()
        if isinstance(extraction_table, pl.DataFrame):
            extraction_table = extraction_table.to_pandas()
        gff_table_bed = BedTool.from_dataframe(gff_table)
        extraction_table_bed = BedTool.from_dataframe(extraction_table)

        if faidx is not None:
            gff_table_bed = gff_table_bed.sort(faidx=faidx)
            extraction_table_bed = extraction_table_bed.sort(faidx=faidx)
        else:
            gff_table_bed = gff_table_bed.sort()
            extraction_table_bed = extraction_table_bed.sort()

        coverage_df = pl.read_csv(
                                gff_table_bed.coverage(extraction_table_bed).fn,
                                has_header=False,
                                separator="\t",
                                new_columns=gff_table.columns.tolist() + GFFExtractor.COVERAGE_FIELDS
                                )\
                        .with_columns(
                                coverage=(1e6 * pl.col("coverage")),
                                atLeastOne=(pl.col("totalHits") > 0).cast(pl.Int32)
                        )\
                        .with_columns(
                                (pl.col("atLeastOne") * pl.col("counts")).alias("atLeastOneUnmerged")
                        )
        if group:
            coverage_df = coverage_df.group_by("compartment",
                                               maintain_order=True)\
                                    .agg(
                                            pl.col("totalHits").sum(),
                                            pl.col("compartmentLength").count().alias("totalCompartments"),
                                            pl.col("counts").sum().alias("totalCompartmentsUnmerged"),
                                            pl.col("overlappingBp").sum(),
                                            pl.col("compartmentLength").sum(),
                                            pl.col("coverage").mean().alias("avg_coverage"),
                                            pl.col("coverage").median().alias("median_coverage"),
                                            pl.col("coverage").std().alias("std_coverage"),
                                            pl.col("coverage").quantile(0.975).alias("q975_coverage"),
                                            pl.col("coverage").quantile(0.025).alias("q025_coverage"),
                                            pl.col("coverage").min().alias("min_coverage"),
                                            pl.col("coverage").max().alias("max_coverage"),
                                            pl.col("atLeastOne").sum().alias("atLeastOne"),
                                            pl.col("atLeastOneUnmerged").sum().alias("atLeastOneUnmerged"),
                                            (1e2 * pl.col("atLeastOne").mean()).alias("perc_atLeastOne"),
                                            (1e2 * pl.col("atLeastOneUnmerged").mean()).alias("perc_atLeastOneUnmerged"),
                                        )\
                                    .with_columns(
                                            (1e6 * pl.col("overlappingBp") / pl.col("compartmentLength")).alias("coverage"),
                                        )
        return coverage_df

    def parse_introns(self, gff_table: pl.DataFrame) -> pl.DataFrame:
        if "compartmentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_compartment_ID(gff_table)
        if "parentID" not in gff_table.columns:
            gff_table = GFFExtractor.parse_parent_ID(gff_table)
        
        gff_table = gff_table.join(
                            gff_table,
                            left_on="parentID",
                            right_on="compartmentID",
                            how="left",
                            )
        return gff_table

    def parse_UTR(self, gff_table: pl.DataFrame) -> pl.DataFrame:
        raise NotImplementedError()

def parse_fasta(fasta: str) -> Iterator[tuple[str, str]]:
    fasta = Path(fasta).resolve()
    if fasta.name.endswith(".gz"):
        f = gzip.open(fasta, "rt")
    else:
        f = open(fasta, mode="r", encoding="utf-8")
    for record in SeqIO.parse(f, "fasta"):
        yield str(record.id), str(record.seq)
    f.close()

def extract_id(accession: str) -> str:
    accession = Path(accession).name
    if accession.count('_') > 2:
        return '_'.join(accession.split('_')[:2])
    return accession.split('.gff')[0]
    # return '.'.join(accession.split('.')[:2])

def extract_name(accession: str, suffix: str = '.tsv') -> str:
    return Path(accession).name.split(f'{suffix}')[0]

class CoverageExtractor:

    def __init__(self, out: str,
                       schedule: str, 
                       design: str, 
                       biotypes: Optional[list[str]] = None,
                       compartments: Optional[list[str]] = None,
                       faidx: Optional[str] = None) -> None:
        self.out = Path(out).resolve()
        self.schedule = Path(schedule).resolve()
        if not self.schedule.is_file():
            raise FileNotFoundError(f"Could not detect schedule file at `{self.schedule}`.")
        self.faidx = faidx

        if compartments is None:
            self.compartments = ["gene", "CDS", "exon", "region"]
        else:
            self.compartments = compartments
        if biotypes is None:
            # protein coding biotype --> parses only protein coding genes coverage
            # non coding biotype --> parses only non coding genges coverage
            # `.` biotype --> refers to generic information: parses genes, exons, CDS, etc. for total coverage
            self.biotypes = ["protein_coding", "non_coding", "."]
        else:
            self.biotypes = biotypes

        self.extractions = dict()
        def _sniff_delimiter(design: str) -> str:
            """Determines if CSV file is tab or comma delimited."""
            with open(design, mode="r", encoding="utf-8") as f:
                for line in f:
                    break
            return "\t" if line.count("\t") > line.count(",") else ","

        with open(design, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=_sniff_delimiter(design))
            for row in reader:
                accession_id = row['accession_id']
                extraction_file = row['extraction']
                self.extractions[accession_id] = extraction_file
        total_extractions = len(self.extractions)
        color = "red" if total_extractions == 0 else "green"
        print(colored(f"Total extractions {total_extractions} have been loaded.", color))
        if total_extractions == 0:
            raise ValueError(f"Could not detected any extractions from design file `{design}`.")

    def load_bucket(self, bucket_id: int) -> dict[str, list[str]]:
        with open(self.schedule, mode="r", encoding="utf-8") as f:
            bucket = json.load(f)[str(bucket_id)]
            print(colored(f"Total {len(bucket)} records have been loaded (bucket {bucket_id}).", "green"))
            return bucket

    @staticmethod
    def _sniff_delimiter(file: str) -> str:
        with open(file, mode="r", encoding="utf-8") as f:
            for line in f:
                break
        return "\t" if line.count("\t") > line.count(",") else ","

    class _TrackProgress:
        def __init__(self, bucket_id: int, 
                           total_records: int,
                           sleeping_time: float) -> None:
            self.track = 0
            self.total_records = total_records
            self.bucket_id = bucket_id
            self.sleeping_time = sleeping_time

        def start(self) -> None:
            while True:
                progress = self.track * 1e2 / self.total_records
                logging.info(f"Current progress for bucket `{self.bucket_id}`: {progress:.2f}.")
                time.sleep(self.sleeping_time)

    def process_bucket(self, bucket_id: int, 
                            partition_col: Optional[str] = None, 
                            group: bool = True,
                            sleeping_time: float = 200) -> None:
        bucket = self.load_bucket(bucket_id=bucket_id)
        logging.info(f"Processing bucket `{bucket_id}`...")
        reader = GFFExtractor(compartments=self.compartments)
        coverage_df = []
        unique_partitions = None
        tracker = CoverageExtractor._TrackProgress(bucket_id=bucket_id,
                                                    total_records=len(bucket),
                                                    sleeping_time=sleeping_time)
        daemon = threading.Thread(target=tracker.start, daemon=True, name="LoggingDaemon")
        daemon.start()
        for gff_file in bucket:
            tracker.track += 1
            accession_id = extract_id(gff_file)
            extraction_filename = self.extractions.get(accession_id)
            if extraction_filename is None:
                logging.info(f"Failed to find extraction file for accession id `{accession_id}`.")
                continue
            delimiter = CoverageExtractor._sniff_delimiter(extraction_filename)
            extraction_table = pl.read_csv(extraction_filename, separator=delimiter)
            extraction_table = extraction_table.rename({col: col[:1].lower() + col[1:] for col in extraction_table.columns})
            if "seqID" not in extraction_table.columns:
                if "chromosome" in extraction_table.columns:
                    extraction_table = extraction_table.rename({
                                                            "chromosome": "seqID",
                                                        })
            if "seqID" not in extraction_table.columns:
                raise KeyError(f"Invalid column for chromosome ID.")

            selection_items = ["seqID", "start", "end"]
            if partition_col is not None:
                selection_items.append(partition_col)
                if partition_col not in extraction_table.columns:
                    raise KeyError(f"Invalid partition column. `{partition_col}` was not found in the dataframe.")
                unique_partitions = set(extraction_table[partition_col])
            extraction_table = extraction_table.select(selection_items)
            gff_table = reader.read_gff(gff_file)
            for biotype in self.biotypes:
                if biotype != ".":
                    gff_table_temp = gff_table.filter(pl.col("biotype") == biotype)
                else:
                    gff_table_temp = gff_table
                gff_table_merged = reader.merge(gff_table=gff_table_temp, faidx=self.faidx)
                coverage_table = reader.parse_coverage(
                                                        gff_table_merged, 
                                                        extraction_table, 
                                                        group=group, 
                                                        faidx=self.faidx)\
                                    .with_columns(
                                            pl.lit(accession_id).alias("#assembly_accession"),
                                            pl.lit(biotype).alias("biotype")
                                        )
                coverage_df.append(coverage_table)
        coverage_df = pl.concat(coverage_df).sort(by=["#assembly_accession", "compartment", "biotype", "coverage"], 
                                                  descending=True)
        coverage_df.write_csv(f"{self.out}/coverage_bucket_{bucket_id}.txt",
                              separator="\t",
                              include_header=True,
                              float_precision=2
                              )
        logging.info(f"Bucket `{bucket_id}` has been processed succesfully.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""""")
    parser.add_argument("schedule", type=str, default="schedule.json")
    parser.add_argument("--out", type=str, default="coverage_out")
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--design", type=str, default="design.csv")
    parser.add_argument("--faidx", type=str, default=None)
    parser.add_argument("--sleeping_time", type=float, default=200)
    parser.add_argument("--group", type=int, default=1)
    parser.add_argument("--partition_col", type=str, default=None)

    args = parser.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(exist_ok=True)
    if not out.is_dir():
        raise ValueError(f'Failure to create destination directory at `{out}`. Is it a nested structure?')

    schedule = args.schedule
    sleeping_time = args.sleeping_time
    group = bool(args.group)
    bucket_id = args.bucket_id
    design = args.design
    faidx = args.faidx

    Path("biologs").mkdir(exist_ok=True)
    logging.basicConfig(
                        level=logging.INFO,
                        filemode="a+",
                        format="%(asctime)s:%(levelname)s:%(message)s",
                        filename=f"biologs/coverage_{bucket_id}.log"
                        )

    extractor = CoverageExtractor(out=out,
                                  schedule=schedule, 
                                  design=design,
                                  faidx=faidx)
    extractor.process_bucket(bucket_id=bucket_id, 
                             sleeping_time=sleeping_time, 
                             group=group, 
                             partition_col=None)
    print(colored(f"Process has been completed succesfully.", "green"))
