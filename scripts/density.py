import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from pybedtools import BedTool
import pybedtools
from mindi.coverage.gff_clean import GFFCleaner
from mindi.coverage.pwm_density import PWMExtractor, strand_evaluators
from mindi.coverage.windows_maker import WindowMaker
from typing import Optional, Callable
from pathlib import Path

GFF_FIELDS = ["seqID", 
          "start", 
          "end", 
          "strand", 
          "chromosome", 
          "motif_start", 
          "motif_end", 
          "sequence", 
          "overlap"]

def merge_and_read(A_df: pd.DataFrame, 
                   c: Optional[list[str]] = None, 
                   o: Optional[list[str]] = None, 
                   delim: str = ";") -> pd.DataFrame:
    if c is None or o is None:
        c, o = [], []
    merged_bed = BedTool.from_dataframe(A_df)\
                        .sort()\
                        .merge(c=c, o=o, delim=delim)
    return pd.read_table(merged_bed.fn, 
                         header=None,
                         names=["seqID", "start", "end"] + c
                        )

def SUM(A_df: pd.DataFrame) -> int:
    return merge_and_read(A_df)\
            .assign(length=lambda ds: ds['end']-ds['start'])['length'].sum()

def query(A_df: pd.DataFrame, 
          B_df: pd.DataFrame, 
          merge_A: bool = False,
          merge_B: bool = False,
          strand_B: bool = False) -> dict:
    # to bed
    A_bed = BedTool.from_dataframe(A_df)
    if "motif_strand" not in B_df:
        B_df = B_df.rename(columns={"strand": "motif_strand"})
    if "strand" in B_df:
        B_df = B_df.drop(columns=["strand"])
    B_bed = BedTool.from_dataframe(B_df)

    # merge if applicable
    if merge_A:
        A_bed = A_bed.merge()
    if merge_B and not strand_B:
        B_bed = B_bed.merge()
    elif merge_B and strand_B:
        B_bed = B_bed.merge(s=True, c="6", o="distinct")

    # sort bed
    A_bed = A_bed.sort()
    B_bed = B_bed.sort()

    # fetch columns
    A_columns = A_df.columns.tolist()
    B_columns = B_df.columns.tolist()[3:]

    # Calculate intersection between target: A_df and query: B_df
    intersect_df = pd.read_table(
            A_bed.intersect(B_bed, wao=True).fn,
            header=None,
            names=A_df.columns.tolist() \
                    + ["chrom", "motif_start", "motif_end"] \
                    + B_columns \
                    + ["overlap"]
        )
    # total queries & targets
    total_queries = A_bed.count()
    total_targets = B_bed.count()
    
    # Estimate:
    # - total number of compartments from A (unique)
    # - the percentage that have at least one overlap
    aggregated_stats = intersect_df.groupby(A_columns, as_index=False)\
                                    .agg(
                                            # total_matches=("overlap", "count"),
                                            match_bp=("overlap", "sum"),
                                        )
    aggregated_stats.loc[:, "atLeastOne"] = (aggregated_stats["match_bp"] > 0).astype(int)
    at_least_one = aggregated_stats.atLeastOne.value_counts()
    try:
        matched_queries = total_queries - int(at_least_one.loc[0])
    except KeyError:
        matched_queries = int(at_least_one.loc[1])

    intersect_df = intersect_df.query("overlap > 0").copy()
    if strand_B and intersect_df.shape[0] > 0:
        intersect_df.loc[:, "non_template"] = (intersect_df.loc[:, "strand"] == intersect_df.loc[:, "motif_strand"]).astype(int)
        
    #   aggregated_stats = intersect_df.groupby(A_columns + ["non_template"], as_index=False)\
    #                                .agg(
    #                                    match_bp=("overlap", "sum")
    #                                    )\
    #                                .pivot(index=A_columns,
    #                                       columns="non_template",
    #                                       values="match_bp"
    #                                       )\
    #                                .fillna(0.0)\
    #                                .melt(id_vars="non_template", value_vars="match_bp")
    #   aggregated_stats.loc[:, "atLeastOne"] = (aggregated_stats["match_bp"] > 0).astype(int)

    # Filter non-overlapping elements
    # intersect_df = intersect_df.query("overlap > 0")

    # Calculate Number of motifs from B mapped to A
    total_motifs = B_df.shape[0]
    # keep compartments from A succesfully mapped to an element of B
    # drop duplicates to ensure that a query from B has not been mapped to more 
    # than one compartment from A; since here we are only interested to estimate 
    # the percentage of motifs from B that map to AT LEAST one compartment.
    # In the previous question, we examined the percentage of compartments from A 
    # that map to AT LEAST ONE from B.
    matched_targets_df = intersect_df\
            .drop_duplicates(subset=['chrom', 'motif_start', 'motif_end'])
    matched_targets = matched_targets_df.shape[0]

    # >> non template proportion
    if strand_B and intersect_df.shape[0] > 0:
        matched_targets_non_template = round(1e2 * intersect_df.non_template.mean(), 2)
        matched_targets_template = round(1e2 - matched_targets_non_template, 2)
    else:
        matched_targets_non_template, matched_targets_template = None, None


    target_perc = round(1e2 * matched_targets / total_targets, 2) if total_targets > 0 else None
    query_perc = round(1e2 * matched_queries / total_queries, 2) if total_queries > 0 else None
    return {
            "intersect": intersect_df,
            "total_targets": total_motifs,
            "matched_targets": matched_targets,
            "matched_targets_non_template": matched_targets_non_template,
            "matched_targets_template": matched_targets_template,
            "target_perc": target_perc,
            "total_queries": total_queries,
            "matched_queries": matched_queries,
            "query_perc": query_perc,
            }

def distance_from_origin(intersect_df: pd.DataFrame, 
                         polynomial: Callable,
                         window_size: int,
                         step: int = 500) -> tuple[pd.DataFrame, pd.Series, np.ndarray]:
    # intersect_df = intersect_df.query("overlap > 0")
    intersect_df["origin"] = intersect_df["start"] + window_size 
    intersect_df["distance"] = np.minimum(
                            np.abs(intersect_df["motif_start"] - intersect_df["origin"]),
                            np.abs(intersect_df["motif_end"] - 1 - intersect_df["origin"])
                    )
    intersect_df["distance_bin"] = intersect_df["distance"].apply(lambda x: x//step+1)
    counts_per_bin = intersect_df.groupby("distance_bin")\
                                 .agg(totalCounts=("seqID", "count"),)
    x_data = counts_per_bin["distance_bin"]
    y_data = counts_per_bin["totalCounts"]
    params, _ = curve_fit(polynomial, x_data, y_data)
    y_pred = polynomial(x_data, *params)
    return counts_per_bin, x_data, y_pred

def bootstrap(intersect_df: pd.DataFrame, 
              window_size: int,
              N: int = 1_000,
              lower_q: float = 0.025,
              upper_q: float = 0.975,
              ) -> tuple[pd.Series, pd.Series, pd.Series]:
    bootstrapped_df = []
    extractor = PWMExtractor()
    for _ in range(N):
        sample_df = intersect_df.sample(frac=1.0, replace=True)
        density = extractor.extract_density(sample_df, 
                                            window_size=window_size,
                                            return_array=True,
                                            enrichment=True
                                            )
        bootstrapped_df.append(density)
    bootstrapped_df = pd.DataFrame(bootstrapped_df)
    average = bootstrapped_df.mean()
    ci_lower = bootstrapped_df.quantile(lower_q)
    ci_upper = bootstrapped_df.quantile(upper_q)
    return average, ci_lower, ci_upper

def relative_density(A_file: str, 
                     B_file: str,
                     window_size: int = 1_000,
                     enrichment: bool = True,
                     ):
    maker = WindowMaker(base=0, window_size=window_size)
    extractor = PWMExtractor()
    A_df = pd.read_table(A_file)
    if "strand" not in A_df:
        A_df["strand"] = "+"
    A_df = A_df[["seqID", "start", "end", "strand"]]
    A_df_win = maker.make_windows(A_df, loci="mid", genome=None)
    B_df = pd.read_table(B_file)
    query_result = query(A_df_win, B_df)
    intersect_df = query_result.pop("intersect")
    density = extractor.extract_density(intersect_df, 
                                        window_size=window_size,
                                        enrichment=enrichment,
                                        return_frame=True
                                        )
    return density

def sniff_delimiter(file: str) -> str:
    with open(file, mode="r", encoding="UTF-8") as f:
        for line in f:
            line = line.strip()
            break
    if "\t" in line:
        return "\t"
    return ","

def determine_strand(seq: str) -> str:
    seq = seq.upper()
    if seq.count("G") >= seq.count("C"):
        return "+"
    return "-"

def extract_density(extraction: str,
                    gff_file: str,
                    window_size: int,
                    mode: str = "density",
                    enrichment: bool = True,
                    compartment: str = "Gene",
                    determine_strand: Optional[Callable[[str], str]] = None,
                    attribute_col: Optional[str] = None,
                    attribute: Optional[str] = "all",
                    tempdir: Optional[str] = None,
                    biotype: Optional[str] = None,
                    genome: Optional[str] = None,
                    ) -> Optional[dict]:
    global GFF_FIELDS
    if tempdir is not None:
        pybedtools.helpers.set_tempdir(tempdir)
        pybedtools.set_tempdir(tempdir)
    if mode != "template" and mode != "density":
        raise ValueError(f"Invalid mode `{mode}` detected.")
    maker = WindowMaker(base=0, window_size=window_size)
    cleaner = GFFCleaner(valid_compartments=[compartment])
    delimiter = sniff_delimiter(extraction)
    try:
        extractions_df = pd.read_csv(extraction, sep=delimiter, comment="#")
        extractions_df.columns = [col[0].lower() + col[1:] for col in extractions_df.columns]
    except pd.errors.EmptyDataError:
        return None
    # change column names to match
    if "seqID" not in extractions_df:
        extractions_df = extractions_df.rename(columns={"chromosome": "seqID", 
                                                        "Chromosome": "seqID", 
                                                        "sequence_name": "seqID", 
                                                        "sequenceName": "seqID"})
    if "seqID" not in extractions_df:
        raise KeyError(f"Column `seqID` is not present in the extractions dataframe ({extraction}).")
    if "end" not in extractions_df:
        extractions_df = extractions_df.rename(columns={"stop": "end", 
                                                        "Stop": "end", 
                                                        "End": "end"})
    if "end" not in extractions_df:
        raise KeyError(f"Column `end` is not present in the extractions dataframe ({extraction}).")
    if "strand" not in extractions_df and "sequence" not in extractions_df and mode == "template":
        raise KeyError(f"Either `strand` or `sequence` must be present in the dataframe when selected mode is `template`.")
    if "strand" not in extractions_df and mode == "template":
        # resolve strand if sequence is present when mode is `template`
        extractions_df["strand"] = extractions_df["sequence"].str.lower().apply(determine_strand)

    if attribute != "all" and attribute_col:
        if attribute_col not in extractions_df:
            raise KeyError(f"Column `{attribute_col}` is not present in the extractions dataframe ({extraction}).")
        extractions_df = extractions_df[extractions_df[attribute_col] == attribute]

    if mode == "template":
        extractions_df = extractions_df[["seqID", "start", "end", "sequence", "strand"]]\
                                .rename(columns={"strand": "motif_strand"})
    else:
        extractions_df = extractions_df[["seqID", "start", "end"]]

    gff_df = cleaner.read_gff(gff_file, 
                              post_filter=[compartment],
                              change_names=True,
                              biotype=biotype is not None, 
                              parse_name=True)
    if biotype:
        gff_df = gff_df.query(f"biotype == '{biotype}'")
        gff_df = gff_df[["seqID", "start", "end", "biotype", "name", "strand"]]
    else:
        gff_df = gff_df[["seqID", "start", "end", "name", "strand"]]
    
    # handle that some GFF files may not have genes or other compartments
    if gff_df.shape[0] == 0:
        # in this case return None, since we cannot make a statistical assesement abou the density
        # there are no compartments to make such assesement
        return None
   
    # return query statistics for template & non template query matches
    if mode == "template":
        strand_B = True
    else:
        strand_B = False

    # Transcription Start Site (TSS) from GFF
    gff_tss = maker.make_windows(gff_df, loci="start", genome=genome)
    tss_query = query(gff_tss, extractions_df, strand_B=strand_B)
    intersect_tss = tss_query.pop("intersect")

    # REPLACE with query function >>>>
    # gff_tss_bed = BedTool.from_dataframe(gff_tss)
    # intersect_df_tss = pd.read_table(
    #                         gff_tss_bed.intersect(extractions_bed, wo=True).fn,
    #                         header=None,
    #                         names=GFF_FIELDS)
    # REPLACE with query function <<<<

    # Transcription End Site (TES) from GFF
    gff_tes = maker.make_windows(gff_df, loci="end", genome=genome)
    tes_query = query(gff_tes, extractions_df, strand_B=strand_B)
    intersect_tes = tes_query.pop("intersect")

    # REPLACE with query function >>>>
    # gff_tes_bed = BedTool.from_dataframe(gff_tes)
    #    intersect_df_tes = pd.read_table(
    #                            gff_tes_bed.intersect(extractions_bed, wo=True).fn,
    #                            header=None,
    #                            names=GFF_FIELDS
    #                            )
    # REPLACE WITH QUERY FUNCTION <<<<

    site_queries = {"tss": tss_query, "tes": tes_query}
    def _extract_density(intersect_df: pd.DataFrame, mode: str = "density") -> pd.DataFrame:
        pwm = PWMExtractor()
        if mode == "template":
            assert determine_strand is not None, f"Determine strand callable was not provided."
            if intersect_df.shape[0] > 0:
                intersect_df["sequence"] = intersect_df["sequence"].str.lower()
                intersect_df["motif_strand"] = intersect_df["sequence"].apply(determine_strand)
            density = pwm.extract_template_density(intersect_df, 
                                                   window_size=window_size, 
                                                   enrichment=enrichment,
                                                )
        elif mode == "density":
            density = pwm.extract_density(intersect_df, 
                                          window_size=window_size, 
                                          return_frame=True,
                                          enrichment=enrichment)
        return density

    densities_df = []
    stats_df = []
    for site, intersect_df in zip(["tss", "tes"], [intersect_tss, intersect_tes]):
        density = _extract_density(intersect_df, mode)
        density = density.astype(int)
        density = density.set_index("index").T
        density.loc[:, "site"] = site
        densities_df.append(density)
        queries_df = pd.DataFrame(site_queries[site], index=[0]).T\
                                    .rename(columns={0: site}).T
        queries_df["matched_queries"] = queries_df["matched_queries"].astype(int)
        queries_df["matched_targets"] = queries_df["matched_targets"].astype(int)
        queries_df["total_queries"] = queries_df["total_queries"].astype(int)
        queries_df["total_targets"] = queries_df["total_targets"].astype(int)
        stats_df.append(queries_df)

    densities_df = pd.concat(densities_df)
    stats_df = pd.concat(stats_df, axis=0)
    return densities_df, stats_df


if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        window_size = int(sys.argv[1])
    else:
        window_size = 500
    density_type = "density"
    determine_strand = strand_evaluators["HDNA"]
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
    parent = Path(__file__).parent.parent
    extracted_file = [file for file in parent.joinpath("extractions").glob("*.tsv") if "IR" in file.name and "GCF_004355385.1" in file.name][0]
    print(extracted_file)
    gff_files = {file.name.split('.gff')[0]: file for file in parent.joinpath("accessions").glob("*.gff")}
    chosen_file_id = extract_id(extracted_file)
    associated_gff = gff_files[chosen_file_id]
    densities = dict()
    # if the name is the same it will overwrite the class above watch out?
    for eDensity in extract_density(
            extraction=extracted_file,
            gff_file=associated_gff,
            window_size=window_size,
            density_type=density_type,
                              determine_strand=determine_strand,
                              ):
        densities.update({(eDensity.loci, eDensity.category): eDensity.density}) 
