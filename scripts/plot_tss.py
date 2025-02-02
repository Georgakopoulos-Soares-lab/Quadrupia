if __name__ == "__main__":

    from plot_enrichment import DensityPlotter, to_png
    from pathlib import Path
    from tqdm import tqdm
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument("--tss", type=str, default="enrichment_compartments.TSS.density.parquet")
    parser.add_argument("--groupby", type=int, default=1)

    args = parser.parse_args()
    tss = args.tss
    groupby = args.groupby
    tss_figures = Path("tss_figures")
    tss_figures.mkdir(exist_ok=True)

    df_tss = pd.read_parquet(tss).drop(columns=['biotype', 'generic'])
    if groupby:
        df_tss = df_tss.groupby(df_tss.index).sum()

    plotter = DensityPlotter()
    index_tss = df_tss.index 
    for i in tqdm(range(df_tss.shape[0])):
        density_tss = df_tss.iloc[i]
        density_tss = density_tss / density_tss.mean()
        index = index_tss[i]
        if density_tss.sum() == 0:
            print(f"Skipping {index}")
            continue

        plotter.plot(density=density_tss)
        filename = f"{tss_figures}/{index}.tss.density.png"
        to_png(plotter.p, filename=filename)
