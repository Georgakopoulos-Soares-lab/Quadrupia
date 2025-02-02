if __name__ == "__main__":

    from plot_enrichment import DensityPlotter, to_png
    from pathlib import Path
    # from bokeh.io import show, output_file
    from tqdm import tqdm
    from bokeh.layouts import column
    import argparse
    import pandas as pd
    import json
    # import polars as pl

    parser = argparse.ArgumentParser(description="""""")
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

    def plot_domain(tss: str, tes: str, domain: str, color: str):
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
            assert df_tss_temp.shape[0] == 3
            tss_average, tss_ci_lower, tss_ci_upper = _unpack(df_tss_temp)

            # tes 
            assert df_tes_temp.shape[0] == 3
            tes_average, tes_ci_lower, tes_ci_upper = _unpack(df_tes_temp)

            plotter.plot_tss_tes(
                         density_tss=tss_average,
                         density_tes=tes_average,
                         ci_lower_tss=tss_ci_lower,
                         ci_upper_tss=tss_ci_upper,
                         ci_lower_tes=tes_ci_lower,
                         ci_upper_tes=tes_ci_upper,
                         use_title=True
                         )
            plots[biotype] = plotter.p
        return plots
    
    figures = {}
    biotypes = []
    def load_colors() -> dict:
        with open("taxonomic_colors.json", mode="r") as f:
            colors = json.load(f)
        return colors

    colors = load_colors()
    # colors = colors["template"]
    colors = colors["non-template"]

    for group in groups:
        tss = tss_files[group]
        tes = tes_files[group]
        plots = plot_domain(tss, tes, group, 
                            color=colors[group] if palette is None else colors[palette],
                            )
        for biotype, p in plots.items():
            figures[group, biotype] = p
            biotypes.append(biotype)
   
    if plot_domains:
        figures_final = {}
        for biotype in biotypes:
            figure = column(figures["Eukaryota", biotype],
                            figures["Archaea", biotype],
                            figures["Viruses", biotype],
                            figures["Bacteria", biotype]
                            )
            figures_final[biotype] = figure 
        to_png(figures_final["proteinCoding"], filename=f"{outdir}/domain_proteinCoding_density.density.png")
        to_png(figures_final["nonCoding"], filename=f"{outdir}/domain_nonCoding_density.density.png")
    else:
        for biotype in tqdm(biotypes, leave=True):
            for group in tqdm(groups, leave=True, position=0):
                to_png(figures[group, biotype], filename=f"{outdir}/{group}_{biotype}_density.density.png")
