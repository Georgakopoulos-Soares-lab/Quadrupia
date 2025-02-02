from termcolor import colored
import csv
from bootstrap_enrichment import Bootstrapper

out = Path(config['out']).resolve()
out.mkdir(exist_ok=True)
TOTAL_BUCKETS = int(config['buckets'])
mode = config['mode']
alpha = round(float(config['alpha']), 2)
DESIGN = config['DESIGN']

# load phylums
PHYLUMS = set()
with open(DESIGN, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter=",")
    for row in reader:
        PHYLUMS.add(row['phylum'].replace(' ', '-'))
PHYLUMS = list(filter(len, PHYLUMS))
total_phylums = len(PHYLUMS)
info_color = "red" if total_phylums == 0 else "blue"
print(colored(f"Total phylums detected: {total_phylums}.", info_color))
# <<

# constants
BIOTYPES = config['biotypes']
if not BIOTYPES:
    BIOTYPES = ["protein_coding", "non_coding"]
DOMAINS = ["Archaea", "Bacteria", "Eukaryota", "Viruses"]
# DOMAINS = ["Bacteria"]
SITES = ["TSS", "TES"]
# <<

# create directories
dest_dir_domain = Path(f"{out}/{mode}/enrichment/domain")
dest_dir_domain.mkdir(exist_ok=True)
dest_dir_phylum = Path(f"{out}/{mode}/enrichment/phylum")
dest_dir_phylum.mkdir(exist_ok=True)

print(f"CHOSEN MODE: `{mode}`.")
print(f"Biotypes: `{BIOTYPES}`.")
print(f"Redirecting domain level outputs to --> `{dest_dir_domain}`.")
print(f"Redirecting phylum level outputs to --> `{dest_dir_phylum}`.")
# <<


rule all:
    input:
        expand(['%s/%s/enrichment/domain/enrichment_bootstrap_alpha_%s.{site}.%s.domain.{domain}.csv' % (out, mode, alpha, mode),
               '%s/%s/enrichment/domain/enrichment_phylums.{site}.%s.{domain}.csv' % (out, mode, mode)], 
                                 site=SITES, domain=DOMAINS),
        expand(['%s/%s/enrichment/phylum/enrichment_bootstrap_alpha_%s.{site}.%s.phylum.{phylum}.csv' % (out, mode, alpha, mode)], site=SITES, phylum=PHYLUMS),

rule taxonomyDomainBootstrap:
    input:
        DESIGN,
        '%s/%s/enrichment/enrichment_compartments.{site}.%s.parquet' % (out, mode, mode),
    output:
        '%s/%s/enrichment/domain/enrichment_bootstrap_alpha_%s.{site}.%s.domain.{domain}.csv' % (out, mode, alpha, mode),
        '%s/%s/enrichment/domain/enrichment_phylums.{site}.%s.{domain}.csv' % (out, mode, mode)
    params:
        window_size=int(config['window_size']),
        alpha=round(float(config['alpha']), 2),
        mode=config['mode'],
        N=int(config['N']),
        join_templates=config['join_templates']
    run:
        bootstrapper = Bootstrapper(design=input[0], enrichment_file=input[1], params=params)
        if params.mode == "template":
            bootstrapper.bootstrap_enrichment(taxonomic_rank="superkingdom", 
                                              rank=wildcards.domain, 
                                              output=output[0])
        else:
            bootstrapper.bootstrap_enrichment_density(taxonomic_rank="superkingdom", 
                                                      rank=wildcards.domain, 
                                                      output=output[0])
        bootstrapper.average_phylums(domain=wildcards.domain, 
                                     output=output[1],
                                     join_templates=params.join_templates)

rule taxonomyPhylumBootstrap:
    input:
        DESIGN,
        '%s/%s/enrichment/enrichment_compartments.{site}.%s.parquet' % (out, mode, mode),
    output:
        '%s/%s/enrichment/phylum/enrichment_bootstrap_alpha_%s.{site}.%s.phylum.{phylum}.csv' % (out, mode, alpha, mode),
    params:
        window_size=int(config['window_size']),
        alpha=round(float(config['alpha']), 2),
        N=int(config['N']),
        mode=config['mode'],
    run:
        bootstrapper = Bootstrapper(design=input[0], enrichment_file=input[1], params=params)
        if params.mode == "template":
            bootstrapper.bootstrap_enrichment(taxonomic_rank="phylum", 
                                              rank=wildcards.phylum.replace('-', ' '), 
                                              output=output[0])
        else:
            bootstrapper.bootstrap_enrichment_density(taxonomic_rank="phylum", 
                                                      rank=wildcards.phylum.replace('-', ' '), 
                                                      output=output[0])
