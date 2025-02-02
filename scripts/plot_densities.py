import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly import graph_objs
import polars as pl
from termcolor import colored

class DensityPlotter(object):
    
    def __init__(self, 
                 tss_path: str, 
                 tes_path: str, 
                 window_size: int = 500,
                 color: str = "lightgreen",
                 lw: float = 1.7,
                 width: int = 1600,
                 height: int = 500) -> None:
        self.window_size = window_size 
        self.tss_path = tss_path 
        self.tes_path = tes_path 
        self.xrange = range(-self.window_size, self.window_size+1)
        self.lw = lw
        self.width = width 
        self.height = height 
        self.color = color
        self.tss_densities: pl.DataFrame 
        self.tes_densities: pl.DataFrame
        self._load_densities()
        print(colored("Dataframes have been loaded.", "green"))

    def _load_densities(self):
        self.tss_densities = pl.read_parquet(self.tss_path).drop(['biotype', 'generic'])
        self.tes_densities = pl.read_parquet(self.tes_path).drop(['biotype', 'generic'])

    def get_tss_tes(self, assembly_id: str) -> tuple:
        df_tss = self.tss_densities.filter(pl.col("#assembly_accession") == assembly_id).drop(['#assembly_accession']).sum()
        df_tes = self.tes_densities.filter(pl.col("#assembly_accession") == assembly_id).drop(['#assembly_accession']).sum()
        df_tss = df_tss / df_tss.mean_horizontal()
        df_tes = df_tes / df_tes.mean_horizontal()
        return list(df_tss.row(0)), list(df_tes.row(0))


    def plot(self, assembly_id: str) -> graph_objs._figure.Figure:
        fig_tss = go.Figure()
        fig_tes = go.Figure()
        df_tss, df_tes = self.get_tss_tes(assembly_id)
        fig_tss.add_trace(
                        go.Scatter(
                            x=list(self.xrange), 
                            y=df_tss,
                            mode='lines', 
                            name="TSS",
                            line=dict(color=self.color, width=self.lw),
                        )
                )
        fig_tes.add_trace(
                    go.Scatter(
                            x=list(self.xrange), 
                            y=df_tes,
                            name="TES",
                            mode='lines', 
                            line=dict(color=self.color, width=self.lw)
                         )
                    )
        fig = make_subplots(
            rows=1, 
            cols=2, 
            subplot_titles=("Transcription Start Site", "Transcription End Site"),
            horizontal_spacing=0.08 
        )
        fig.add_trace(fig_tss.data[0], row=1, col=1)
        fig.add_trace(fig_tes.data[0], row=1, col=2)
        ymax_tss = max(fig_tss.data[0].y)
        ymax_tes = max(fig_tes.data[0].y)
        fig.add_shape(
                    type='line', 
                    x0=0, 
                    y0=0, 
                    x1=0, 
                    y1=ymax_tss, 
                    line=dict(color='black', width=3, dash='dash'),
                    layer='above',
                    xref="x1", 
                    yref="y1", 
                    row=1, 
                    col=1)
        fig.add_shape(
                    type='line', 
                    x0=-self.window_size, 
                    y0=1.0, 
                    x1=self.window_size, 
                    y1=1.0, 
                    line=dict(color='crimson', width=2), #, dash='dash'),
                    # layer='above',
                    xref="x1", 
                    yref="y1", 
                    row=1, 
                    col=1)
        fig.add_shape(
                    type='line', 
                    x0=0, 
                    y0=0, 
                    x1=0, 
                    y1=ymax_tes, 
                    xref="x2", 
                    yref="y2", 
                    row=1, 
                    col=2,
                    line=dict(color='black', width=3, dash='dash'),
                    layer='above',
                    )
        fig.add_shape(type='line', 
                    x0=-self.window_size, 
                    y0=1.0, 
                    x1=self.window_size, 
                    y1=1.0, 
                    line=dict(color='crimson', width=2),# dash='dash'),
                    # layer='above',
                    xref="x2", 
                    yref="y2", 
                    row=1, 
                    col=2)
        fig.update_xaxes(
            tickmode='linear', 
            tick0=-500, 
            dtick=100, 
            title='',
            row=1, 
            tickfont=dict(size=16),
            col=1
        )
        fig.update_yaxes(
            title="Enrichment", 
            range=[0, None], 
            title_font=dict(size=24),
            tickfont=dict(size=16),
            row=1, 
            col=1
        )
        fig.update_xaxes(
            tickmode='linear', 
            tick0=-500, 
            dtick=100, 
            tickfont=dict(size=16), 
            title='', 
            row=1,
            col=2
        )
        fig.update_yaxes(
            title="Enrichment", 
            title_font=dict(size=24),
            range=[0, None], 
            tickfont=dict(size=18),
            row=1, 
            col=2
        )
        fig.update_layout(
            title="",
            title_x=0.5,
            # title_font=dict(family="Arial", size=34, color="darkblue"),
            title_xanchor="center",
            plot_bgcolor="white", 
            paper_bgcolor="whitesmoke",
            hovermode="closest",
            yaxis=dict(
                showgrid=True,    
                range=[0, None],
                gridcolor='lightgray', 
                gridwidth=1,          
                tick0=0,  
                dtick=0.5,             
                title_font=dict(size=24),
                tickfont=dict(size=16)   
            ),
            yaxis2=dict(
                showgrid=True,          
                gridcolor='lightgray', 
                range=[0, None],
                gridwidth=1,             
                tick0=0,              
                dtick=0.5,               
                title_font=dict(size=24),
                tickfont=dict(size=16)    
            ),
            annotations=[
                    dict(
                        x=0.23, 
                        y=1.02, 
                        xref="paper", 
                        yref="paper", 
                        text="Transcription Start Site", 
                        showarrow=False,
                        font=dict(size=26) 
                    ),
                    dict(
                        x=0.77, 
                        y=1.02, 
                        xref="paper", 
                        yref="paper", 
                        text="Transcription End Site", 
                        showarrow=False,
                        font=dict(size=26)
                    )
            ],
            showlegend=False,
            margin=dict(l=40, r=40, t=40, b=40),
            width=self.width,
            height=self.height, 
        )
        return fig


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="""Specify parameters sto customize the figure.""")
    parser.add_argument("--assembly_id", type=str, default="GCF_000181295.1")
    parser.add_argument("--tss_path", type=str, default="../../enrichment_compartments.TSS.density.parquet")
    parser.add_argument("--tes_path", type=str, default="../../enrichment_compartments.TES.density.parquet")
    parser.add_argument("--color", type=str, default="#8DD3C7")
    parser.add_argument("--lw", type=float, default=2.3)
    parser.add_argument("--width", type=int, default=1500)
    parser.add_argument("--height", type=int, default=500)
    args = parser.parse_args()

    assembly_id = args.assembly_id
    color = args.color 
    width = args.width 
    height = args.height 
    lw = args.lw 
    tss_path = args.tss_path
    tes_path = args.tes_path

    # plot
    density = DensityPlotter(tss_path=tss_path, 
                             tes_path=tes_path,
                             color=color,
                             height=height,
                             width=width,
                             lw=lw)
    fig = density.plot(assembly_id=assembly_id)
    fig.show()
