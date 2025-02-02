if __name__ == "__main__":

    from plot_enrichment import DensityPlotter, to_png
    from pathlib import Path
    # from bokeh.io import show, output_file
    from tqdm import tqdm
    from bokeh.layouts import column
    import argparse
    import pandas as pd
    import json
    import polars as pl

    parser = argparse.ArgumentParser(description="""Plots figures for TSS/TES Densities.""")
    parser.add_argument("--indir", type=str, default="domain")
    parser.add_argument("--outdir", type=str, default="tss_tes_figures")
    parser.add_argument("--plot_domains", type=int, default=1)
    parser.add_argument("--palette", type=str)
    parser.add_argument("--rank", type=str, default="superkingdom")

    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)
    plot_domains = args.plot_domains
    palette = args.palette
    rank = args.rank

    tree = pl.read_csv("tree_of_life_new.txt.gz", separator="\t")\
                            .with_columns(pl.col(rank).str.replace_all(" ", "-", literal=True))
    tree_mapper = dict(zip(tree[rank], tree["superkingdom"]))
    def load_colors() -> dict:
        with open("taxonomic_colors.json", mode="r") as f:
            colors = json.load(f)
        return colors
    colors = load_colors()

    domains = ["Eukaryota", "Archaea", "Bacteria", "Viruses"]
    density_files = [file for file in indir.glob("*.csv")]
    tss_files = {file.name.split('.')[-2]: file for file in density_files if "TSS" in file.name and file.name.startswith("enrichment_bootstrap")}
    tes_files = {file.name.split('.')[-2]: file for file in density_files if "TES" in file.name and file.name.startswith("enrichment_bootstrap")}
    groups = tss_files.keys()

    def load_density(site: str) -> pd.DataFrame:
        df = pd.read_csv(site, index_col=0)
        df["typ"] = df.index.map(lambda x: x.split("_")[2])
        df["biotype"] = df.index.map(lambda x: x.split("_")[3])
        return df

    def plot_domain(tss: str, tes: str, domain: str, color: str, template_color: str):
        def _unpack(df: pd.DataFrame) -> tuple:
            return df.iloc[0], df.iloc[1], df.iloc[2]
        df_tss = load_density(tss)
        df_tes = load_density(tes)
        plotter = DensityPlotter(title=domain, 
                                 color=color, 
                                 color_template=template_color,
                                 confidence_color=color,
                                 confidence_color_template=template_color,
                                 confidence_alpha=0.3,
                                 lw=2.5,
                                 )
        plots = {}
        for biotype in ["proteinCoding", "nonCoding"]:
            df_tss_temp = df_tss[df_tss["biotype"] == biotype].drop(columns=[rank, 'biotype'])
            df_tes_temp = df_tes[df_tes["biotype"] == biotype].drop(columns=[rank, 'biotype'])

            # tes
            df_tss_template = df_tss_temp[df_tss_temp["typ"] == "template"]\
                                        .drop(columns=['typ'])
            df_tss_nontemplate = df_tss_temp[df_tss_temp["typ"] == "non-template"]\
                                        .drop(columns=['typ'])
            assert df_tss_template.shape[0] == 3
            assert df_tss_nontemplate.shape[0] == 3
            tss_average_template, tss_ci_lower_template, tss_ci_upper_template = _unpack(df_tss_template)
            tss_average_nontemplate, tss_ci_lower_nontemplate, tss_ci_upper_nontemplate = _unpack(df_tss_nontemplate)

            # tes 
            df_tes_template = df_tes_temp[df_tes_temp["typ"] == "template"]\
                                        .drop(columns=['typ'])
            df_tes_nontemplate = df_tes_temp[df_tes_temp["typ"] == "non-template"]\
                                        .drop(columns=['typ'])
            assert df_tes_template.shape[0] == 3
            assert df_tes_nontemplate.shape[0] == 3
            tes_average_template, tes_ci_lower_template, tes_ci_upper_template = _unpack(df_tes_template)
            tes_average_nontemplate, tes_ci_lower_nontemplate, tes_ci_upper_nontemplate = _unpack(df_tes_nontemplate)

            plotter.plot_tss_tes(
                         density_tss=tss_average_nontemplate,
                         density_tes=tes_average_nontemplate,
                         ci_lower_tss=tss_ci_lower_nontemplate,
                         ci_upper_tss=tss_ci_upper_nontemplate,
                         ci_lower_tes=tes_ci_lower_nontemplate,
                         ci_upper_tes=tes_ci_upper_nontemplate,
                         density_tss_template=tss_average_template,
                         density_tes_template=tes_average_template,
                         ci_lower_tss_template=tss_ci_lower_template,
                         ci_upper_tss_template=tss_ci_upper_template,
                         ci_lower_tes_template=tes_ci_lower_template,
                         ci_upper_tes_template=tes_ci_upper_template,
                         use_title=True
                         )
            plots[biotype] = plotter.p
        return plots
    
    figures = {}
    biotypes = set()
    template_colors = colors["template"]
    nontemplate_colors = colors["non-template"]

    for group in groups:
        tss = tss_files[group]
        tes = tes_files[group]
        palette = tree_mapper[group]
        plots = plot_domain(tss, 
                            tes, 
                            group, 
                            color=nontemplate_colors[group] if palette is None else nontemplate_colors[palette],
                            template_color=template_colors[group] if palette is None else template_colors[palette]
                            )
        for biotype, p in plots.items():
            figures[group, biotype] = p
            biotypes.add(biotype)
   
    if plot_domains:
        figures_final = {}
        for biotype in biotypes:
            figure = column(figures["Eukaryota", biotype],
                            figures["Archaea", biotype],
                            figures["Viruses", biotype],
                            figures["Bacteria", biotype]
                            )
            figures_final[biotype] = figure 
        to_png(figures_final["proteinCoding"], filename=f"{outdir}/domain_proteinCoding_density.template.png")
        to_png(figures_final["nonCoding"], filename=f"{outdir}/domain_nonCoding_density.template.png")
    else:
        for biotype in tqdm(biotypes, leave=True):
            for group in tqdm(groups, leave=True, position=0):
                to_png(figures[group, biotype], filename=f"{outdir}/{group}_{biotype}_density.template.png")
