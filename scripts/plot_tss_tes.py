if __name__ == "__main__":

    from plot_enrichment import DensityPlotter, to_png
    from pathlib import Path
    from tqdm import tqdm
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument("--tss", type=str, default="enrichment_compartments.TSS.density.parquet")
    parser.add_argument("--tes", type=str, default="enrichment_compartments.TES.density.parquet")
    parser.add_argument("--groupby", type=int, default=1)

    args = parser.parse_args()
    tss = args.tss
    tes = args.tes
    groupby = args.groupby
    tss_figures = Path("tss_figures")
    tss_figures.mkdir(exist_ok=True)

    df_tss = pd.read_parquet(tss).drop(columns=['biotype', 'generic'])
    if groupby:
        df_tss = df_tss.groupby(df_tss.index).sum()
    df_tes = pd.read_parquet(tes).drop(columns=['biotype', 'generic'])
    if groupby:
        df_tes = df_tes.groupby(df_tes.index).sum()

    plotter = DensityPlotter()
    index_tss = df_tss.index 
    for i in tqdm(range(df_tss.shape[0]), leave=True):
        density_tss = df_tss.iloc[i]
        density_tss = density_tss / density_tss.mean()
        density_tes = df_tes.iloc[i]
        density_tes = density_tes / density_tes.mean()
        index = index_tss[i]
        if density_tss.sum() == 0:
            print(f"Skipping {index}.")
            continue

        plotter.plot_tss_tes(density_tss=density_tss,
                            density_tes=density_tes)
        plotter.show()
        filename = f"{tss_figures}/{index}.tss-tes.density.png"
        to_png(plotter.p, filename=filename)
