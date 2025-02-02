import pandas as pd 
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass, field
import json
from termcolor import colored
from tqdm import tqdm
from typing import Optional
from scipy.stats import ks_2samp

from mindi.coverage.pwm_density import PWMExtractor

# plotting
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass(frozen=True)
class params:
    N: int = field(default=1000)
    window_size: int = field(default=500)
    alpha: float = field(default=0.05)

def bootstrap_density(intersect_df: pd.DataFrame, 
                        window_size: int, 
                        nsamples: int = 1000, 
                        alpha: float = 0.05) -> tuple:
    df_bootstrap = []
    extractor = PWMExtractor()
    for _ in tqdm(range(nsamples), leave=True, position=0):
        df_sample = intersect_df.sample(frac=1.0, replace=True)
        df_sample = extractor.extract_density(df_sample, 
                                              window_size=window_size, 
                                              enrichment=True,
                                              return_array=True)
        df_bootstrap.append(df_sample)
    df_bootstrap = pd.DataFrame(df_bootstrap)
    mean = df_bootstrap.mean()
    ci_lower = df_bootstrap.quantile(alpha/2)
    ci_upper = df_bootstrap.quantile(1-alpha/2)
    return mean, ci_lower, ci_upper

def bootstrap(df: pd.DataFrame, 
              N: int = 1_000, 
              alpha: float = 0.05) -> tuple:
    bootstrapped_samples = []
    for _ in tqdm(range(N), leave=True, position=0):
        sample_df = df.sample(frac=1.0, replace=True)
        sample_df = sample_df.sum(axis=0)
        mean = sample_df.mean()
        resampled_sum = sample_df / mean
        bootstrapped_samples.append(resampled_sum)
    bootstrapped_samples = pd.DataFrame(bootstrapped_samples)
    average = bootstrapped_samples.mean(axis=0)
    # two-tailed interval (1-a)%
    lower_ci = bootstrapped_samples.quantile(alpha/2)
    upper_ci = bootstrapped_samples.quantile(1 - alpha/2)
    return average, lower_ci, upper_ci

class Bootstrapper:
    
    def __init__(self, enrichment_file: str, design: str, params=params) -> None:
        self.enrichment_df = None
        self.enrichment_file = Path(enrichment_file).resolve()
        self.design = Path(design).resolve()
        self.params = params
        self.combinations = [
                        ("Occurrences", "protein_coding"),
                        ("Occurrences", "non_coding") ]
        self.taxonomic_ranks = ["phylum", "kingdom", "superkingdom"]
        if not self.design.is_file():
            raise FileNotFoundError(f"Could not detect design file `{design}`.") 
        if not self.enrichment_file.is_file():
            raise FileNotFoundError(f"Could not detect enrichment file `{enrichment_file}`.") 

    def load_table(self, taxonomic_ranks: Optional[list[str]] = None):
        if self.enrichment_df is not None:
            return self
        if taxonomic_ranks is None:
            taxonomic_ranks = self.taxonomic_ranks
        if not isinstance(taxonomic_ranks, list):
            raise TypeError(f"Invalid type for taxonomic ranks. Expected list, but received {type(taxonomic_ranks)}.")
        design_df = pd.read_csv(self.design, usecols=["accession_id"] + taxonomic_ranks)
        enrichment_df = pd.read_parquet(self.enrichment_file, engine="fastparquet")\
                        .merge(
                                design_df,
                                right_on="accession_id",
                                left_on="#assembly_accession",
                                how="inner"
                              )
        self.enrichment_df = enrichment_df
        return self

    def bootstrap_enrichment_density(self, 
                                     taxonomic_rank: str, 
                                     rank: str, 
                                     output: str, 
                                     taxonomic_ranks: Optional[list[str]] = None, 
                                     combinations: Optional[list[tuple]] = None) -> None:

        enrichment_df = self.load_table(taxonomic_ranks=taxonomic_ranks).enrichment_df
        if taxonomic_rank not in enrichment_df:
            raise KeyError(f"Invalid specified taxonomic rank `{taxonomic_rank}`.")
        enrichment_df = enrichment_df.query(f"{taxonomic_rank} == '{rank}'").copy()

        print(f"Initializing bootstrap for taxonomic rank {taxonomic_rank} with value=`{rank}`.")
        print(f"Specified confidence: {1e2 * (1-self.params.alpha):.2f}")
        print(f"Total resampling iterations: {self.params.N}")
        if enrichment_df.shape[0] == 0:
            raise ValueError(f"Empty dataframe for taxonomic rank `{taxonomic_rank}` with value=`{rank}`.")
        if combinations is None:
            combinations = self.combinations
        ## Domain level bootstrap
        confidence_intervals = defaultdict(list)
        for comb in combinations: 
            _, biotype = comb
            temp_df = enrichment_df[enrichment_df["biotype"] == biotype]\
                        [[str(i) for i in range(-params.window_size, params.window_size+1)]]
            average, lower_ci, upper_ci = bootstrap(temp_df, N=params.N, alpha=params.alpha)
            biotype = biotype.replace("_coding", "Coding")
            confidence_intervals[f"average_{biotype}_{rank}"] = average
            confidence_intervals[f"lowerCI_{biotype}_{rank}"] = lower_ci
            confidence_intervals[f"upperCI_{biotype}_{rank}"] = upper_ci
        confidence_intervals = pd.DataFrame(confidence_intervals).T
        confidence_intervals["biotype"] = confidence_intervals.index.map(lambda x: x.split("_")[2])
        confidence_intervals[taxonomic_rank] = rank
        confidence_intervals = confidence_intervals[[taxonomic_rank, "biotype"] + [str(i) for i in range(-params.window_size, params.window_size+1)]]
        confidence_intervals.to_csv(output, sep=",", mode="w", header=True, index=True)
        print(colored(f"Bootstrap has succesfully been completed for taxonomic rank {taxonomic_rank}=`{rank}`.", "green"))
        return
        
    def bootstrap_enrichment(self, taxonomic_rank: str, rank: str, output: str) -> None:
        enrichment_df = self.load_table().enrichment_df
        if taxonomic_rank not in enrichment_df:
            raise KeyError(f"Invalid specified taxonomic rank `{taxonomic_rank}`.")
        enrichment_df = enrichment_df.query(f"{taxonomic_rank} == '{rank}'").copy()

        print(f"Initializing bootstrap for taxonomic rank {taxonomic_rank} with value=`{rank}`.")
        print(f"Specified confidence: {1e2 * (1-self.params.alpha):.2f}")
        print(f"Total resampling iterations: {self.params.N}")
        if enrichment_df.shape[0] == 0:
            raise ValueError(f"Empty dataframe for taxonomic rank `{taxonomic_rank}` with value=`{rank}`.")
        
        combinations = [("Occurrences_template", "protein_coding"),
                        ("Occurrences_template", "non_coding"),
                        ("Occurrences_non_template", "protein_coding"),
                        ("Occurrences_non_template", "non_coding")]
        ## Domain level bootstrap
        confidence_intervals = defaultdict(list)
        for comb in combinations: 
            typ, biotype = comb
            temp_df = enrichment_df[(enrichment_df["template|non_template"] == typ) & (enrichment_df["biotype"] == biotype)]\
                        [[str(i) for i in range(-params.window_size, params.window_size+1)]]
            average, lower_ci, upper_ci = bootstrap(temp_df, N=params.N, alpha=params.alpha)
            biotype = biotype.replace("_coding", "Coding")
            typ = typ.replace("non_template", "non-template").replace("Occurrences_", "")
            confidence_intervals[f"average_{typ}_{biotype}_{rank}"] = average
            confidence_intervals[f"lowerCI_{typ}_{biotype}_{rank}"] = lower_ci
            confidence_intervals[f"upperCI_{typ}_{biotype}_{rank}"] = upper_ci

        #        for biotype in ["protein_coding", "non_coding"]:
        #            avg_template = confidence_intervals[f"average_Occurrences_template_{biotype}_{domain}"]
        #            avg_non_template = confidence_intervals[f"average_Occurrences_non_template_{biotype}_{domain}"]
        #            significance = ks_2samp(avg_template, avg_non_template, alternative='two-sided')
        confidence_intervals = pd.DataFrame(confidence_intervals).T
        confidence_intervals["biotype"] = confidence_intervals.index.map(lambda x: x.split("_")[3])
        confidence_intervals[taxonomic_rank] = rank
        confidence_intervals["typ"] = confidence_intervals.index.map(lambda x: x.split("_")[2])
        confidence_intervals = confidence_intervals[[taxonomic_rank, "biotype", "typ"] + [str(i) for i in range(-params.window_size, params.window_size+1)]]
        confidence_intervals.to_csv(output, sep=",", mode="w", header=True, index=True)
        print(colored(f"Bootstrap has succesfully been completed for taxonomic rank {taxonomic_rank}=`{rank}`.", "green"))
        return

    def average_phylums(self, domain: str, output: str, join_templates: int = 0) -> None:
        self.load_table()
        enrichment_df = self.enrichment_df
        if "superkingdom" not in enrichment_df:
            raise KeyError(f"Invalid specified domain `{domain}`.")
        enrichment_df = enrichment_df.query(f"superkingdom == '{domain}'").copy()
        if "phylum" not in enrichment_df: 
            raise KeyError(f"No phylums detected for domain `{domain}`.")
        ## Phylum calculation
        ## keep assembly accessions with a present phylum
        enrichment_df_phylum = enrichment_df.dropna(subset=['phylum'])
        phylum_to_kingdom = enrichment_df_phylum[["phylum", "kingdom"]].set_index("phylum")["kingdom"].to_dict()
        if "template|non_template" in enrichment_df_phylum and not join_templates:
            print(colored("Template & Non-template partition detected.", "blue"))
            enrichment_df_phylum = enrichment_df_phylum.groupby(["phylum", 
                                                                 "template|non_template", 
                                                                 "biotype"])
        else:
            enrichment_df_phylum = enrichment_df_phylum.groupby(["phylum", "biotype"])
        enrichment_df_phylum = enrichment_df_phylum.agg({str(i): "sum" \
                                                            for i in range(-params.window_size, params.window_size+1)
                                                        })
        enrichment_df_phylum = enrichment_df_phylum.apply(lambda row: row / row.mean(), axis=1)
        for col in range(-params.window_size, params.window_size+1):
            enrichment_df_phylum[str(col)] = enrichment_df_phylum[str(col)].round(2)
        enrichment_df_phylum.reset_index(inplace=True)
        enrichment_df_phylum.loc[:, "kingdom"] = enrichment_df_phylum["phylum"].map(phylum_to_kingdom)
        enrichment_df_phylum.loc[:, "domain"] = domain
        enrichment_df_phylum.to_csv(output, sep=",", mode="w", header=True, index=True)
        print(colored(f"Phylum averaging has succesfully been completed for domain=`{domain}`.", "green"))
        return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Utility""")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--N", type=int, default=1_000)
    parser.add_argument("--design", type=str, default="design.csv")
    parser.add_argument("--enrichment", type=str, default="")
    parser.add_argument("--mode", type=str, default="template", choices=["template", "density"])
    parser.add_argument("--taxonomic_rank", type=str, default="domain", choices=["domain", "kingdom", "phylum"])
    parser.add_argument("--output", type=str, default="bootstrap.txt")
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--rank", type=str, default="Bacteria", choices=["Eukaryota", "Archaea", "Viruses", "Bacteria"])
    args = parser.parse_args()
    N = args.N 
    alpha = args.alpha
    design = args.design 
    rank = args.rank
    window_size = int(args.window_size)
    taxonomic_rank = args.taxonomic_rank
    enrichment_file = args.enrichment
    output = args.output
    mode = args.mode
    param = params(alpha=alpha, N=N, window_size=window_size)
    bootstrapper = Bootstrapper(params=param, 
                                design=design, 
                                enrichment_file=enrichment_file)
    if mode == "template":
        bootstrapper.bootstrap_enrichment(taxonomic_rank=taxonomic_rank,
                                      rank=rank,
                                      output=output)
    else:
        bootstrapper.bootstrap_enrichment_density(taxonomic_rank=taxonomic_rank,
                                                  rank=rank,
                                                  output=output)
