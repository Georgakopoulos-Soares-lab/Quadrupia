from bokeh.io import output_notebook, show, output_file
from bokeh.layouts import gridplot, column, layout, row
from bokeh.plotting import figure
from bokeh.models import Label, Legend, Span, Title, Label, Div, Spacer, LegendItem
import numpy as np
from typing import Optional 
from dataclasses import dataclass
from bokeh.plotting import save 
from bokeh.io.export import export_png
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass(frozen=True, slots=True, kw_only=True)
class LabelModel(object):
    label: str
    size: int 
    tick_size: int

@dataclass(frozen=True, kw_only=True, slots=True)
class PlotObject(object):
    line_width: int
    line_dash: str
    color: str

class DensityPlotter:
    
    def __init__(self, 
                 window: int = 500,
                 width: int = 750, 
                 height: int = 450, 
                 color: str = "black",
                 color_template: str = "lightblue",
                 lw: float = 2.3,
                 opacity: float = 1.0,
                 opacity_template: float = 0.6,
                 xlabel: str = "", 
                 ylabel: str = "Enrichment",
                 xlabel_size: int = 15, 
                 ylabel_size: int = 18,
                 xlabel_tick_size: int = 16,
                 ylabel_tick_size: int = 14,
                 horizontal_line_width: float = 2.3,
                 horizontal_line_dash: str = "dashed",
                 horizontal_color: str = "crimson",
                 vertical_line_width: float = 2.3,
                 vertical_line_dash: str = "dashed",
                 vertical_color: str=  "black",
                 confidence_color: str = "blue",
                 confidence_color_template: str = "lightblue",
                 confidence_alpha: float = 0.2,
                 confidence_alpha_template: float = 0.2,
                 confidence_label: str = "",
                 title: str = "Bokeh Plot",
                 title_font_size: int = 16,
                 title_box_height: int = 50,
                 title_border_width: int = 3,
                 title_font_weight: str = "bold",
                 xgrid_line_dash: str = 'dotted',
                 ygrid_line_dash: str = 'dotted',
                 grid_line_alpha: float = 0.4,
                 background_fill_color: str = '#f0f0f5',
                 border_fill_color: str = '#ffffff',
                 toolbar_location: Optional[str] = None,
                 ) -> None:
        
        # figure attributes
        self.window = window
        self.width = width
        self.height = height 
        self.background_fill_color = background_fill_color
        self.border_fill_color = border_fill_color

        # lineplot
        self.color = color 
        self.color_template = color_template
        self.lw = lw
        self.opacity = opacity 
        self.opacity_template = opacity_template
        self.xlabel = LabelModel(label=xlabel, size=xlabel_size, tick_size=xlabel_tick_size)
        self.ylabel = LabelModel(label=ylabel, size=ylabel_size, tick_size=ylabel_tick_size)

        # confidence interval attributes
        self.confidence_color = confidence_color 
        self.confidence_alpha = confidence_alpha
        self.confidence_label = confidence_label
        self.confidence_color_template = confidence_color_template
        self.confidence_alpha_template = confidence_alpha_template

        # toolbar
        self.toolbar_location = toolbar_location

        # Grid lines attributes
        self.xgrid_line_dash = xgrid_line_dash 
        self.ygrid_line_dash = ygrid_line_dash
        self.grid_line_alpha = grid_line_alpha

        # vertical & horizontal lines
        self.horizontal = PlotObject(line_width=horizontal_line_width,
                                     line_dash=horizontal_line_dash,
                                     color=horizontal_color)
        self.vertical = PlotObject(line_width=vertical_line_width,
                                   line_dash=vertical_line_dash,
                                   color=vertical_color)
        
        # title attributes
        self.title = title
        self.title_font_size = title_font_size
        self.title_box_height = title_box_height
        self.title_border_width = title_border_width
        self.title_font_weight = title_font_weight

        # figure
        self.p = None

    def show(self) -> None:
        show(self.p)

    def plot_tss_tes(self,
                     density_tss,
                     density_tes,
                     ci_lower_tss = None,
                     ci_upper_tss = None,
                     ci_lower_tes = None,
                     ci_upper_tes = None,
                     density_tss_template = None,
                     ci_lower_tss_template = None,
                     ci_upper_tss_template = None,
                     density_tes_template = None,
                     ci_lower_tes_template = None,
                     ci_upper_tes_template = None,
                     use_title: bool = False
                     ):
        p_tss = self.plot(density_tss, 
                          ci_lower=ci_lower_tss,
                          ci_upper=ci_upper_tss,
                          density_template=density_tss_template,
                          ci_lower_template=ci_lower_tss_template,
                          ci_upper_template=ci_upper_tss_template,
                          title="",
                          )
        p_tes = self.plot(density_tes, 
                          ci_lower=ci_lower_tes,
                          ci_upper=ci_upper_tes,
                          density_template=density_tes_template,
                          ylabel="",
                          ci_lower_template=ci_lower_tes_template,
                          ci_upper_template=ci_upper_tes_template,
                          title="",
                          )
        if use_title:
            title_div = Div(
                    text=f"""<div style="width: 100%; height: {self.title_box_height}px; padding: 10px 650px; 
                      font-size: {self.title_font_size}pt; 
                      font-weight: {self.title_font_weight}; 
                      margin-bottom: 0px; 
                      margin-left: 60px; 
                      text-align: center; 
                      display: flex; 
                      justify-content: center;"> """ + f"""{self.title}</div>""",
                    # width=200,
                    # height=60
            )
            p = column(title_div, row(p_tss, Spacer(width=20), p_tes))
        else:
            p = row(p_tss, Spacer(width=20), p_tes)
        self.p = p
        return self.p

    def plot(self, 
             density, 
             ci_lower: Optional[str] = None, 
             ci_upper: Optional[str] = None,
             density_template = None,
             ci_lower_template = None,
             ci_upper_template = None,
             xlabel: Optional[str] = None,
             ylabel: Optional[str] = None,
             title: Optional[str] = None,
             ymin: Optional[float] = None,
             ymax: Optional[float] = None,
             use_title: bool = False,
             ):
        padding: float = 0.5
        # y-axis boundaries
        if ymax is None:
            ymax = max(np.max(density), np.max(ci_upper))
            if density_template is not None:
                ymax = max(np.max(density_template), ymax)
            if ci_upper_template is not None:
                ymax = max(np.max(ci_upper_template), ymax)
            ymax += 0.6
        if ymin is None:
            ymin = min(np.min(density), np.min(ci_lower))
            if density_template is not None:
                ymin = min(np.min(density_template), ymin)
            if ci_lower_template is not None:
                ymin = min(np.min(ci_lower_template), ymin)
            ymin -= 0.6

        p = figure(
                    width=self.width, 
                    height=self.height, 
                    x_range=(-self.window - padding, self.window + padding),
                    y_range=(0, ymax),
                    # min_border=0,
                    toolbar_location=self.toolbar_location
                )
        # indices = list(density.index)
        density.index = density.index.map(int)
        y = density.values 
        x = density.index
        if density_template is None:
            p.line(x, 
               y, 
               line_width=self.lw, 
               color=self.color, 
               alpha=self.opacity,
            ) 
        else:
            x = density.index.map(int)
            y = density.values
            p.line(
                    x,
                    y,
                    line_width=self.lw,
                    color=self.color_template,
                    alpha=self.opacity_template,
                    # legend_label="Template"
                    )
            x = density_template.index.map(int)
            y = density_template.values
            p.line(x, 
                   y, 
                   line_width=self.lw, 
                   color=self.color, 
                   alpha=self.opacity,
                   # legend_label="Non Template"
                ) 
        if ci_lower is not None and ci_upper is not None:
            p.varea(x, 
                    ci_lower, 
                    ci_upper, 
                    fill_color=self.confidence_color, 
                    fill_alpha=self.confidence_alpha, 
                    # legend_label=self.confidence_label
                )
            if ci_lower_template is not None and ci_upper_template is not None:
                p.varea(x, 
                    ci_lower_template, 
                    ci_upper_template, 
                    fill_color=self.confidence_color_template, 
                    fill_alpha=self.confidence_alpha_template, 
                    # legend_label=self.confidence_label
                    )


            rect_glyph = p.rect(x=[0], 
                                y=[0], 
                                width=1, 
                                height=1, 
                                color=self.color_template, 
                                visible=False)
            item_template = LegendItem(label="Template", renderers=[rect_glyph])
            rect_glyph = p.rect(x=[0], 
                                y=[0], 
                                width=1, 
                                height=1, 
                                color=self.color,
                                visible=False)
            item_nontemplate = LegendItem(label="Non Template", renderers=[rect_glyph])
            legend = Legend(items=[item_template, item_nontemplate],  # Add the rectangle legend item
                            title="",
                            title_text_font_size="12pt",
                            label_text_font_size="14pt",
                            glyph_height=20,
                            glyph_width=30,
                            spacing=10)   
            legend.location = "top_right"
            p.add_layout(legend, "center")
                
        # Modify grid lines
        # p.xgrid.grid_line_dash = self.xgrid_line_dash
        # p.ygrid.grid_line_dash = self.ygrid_line_dash
        # p.grid.grid_line_alpha = self.grid_line_alpha

        # p.title = Title(text="Fancy Plot Title", 
        #                 align="center", 
        #                 text_font_size="16pt", 
        #                 text_font_style="bold", 
        #                 text_color="black")
        # p.title_location = "above"
        # p.title.text_font_size = "16pt"
        # p.title.text_font_style = "bold"
        # p.title.text_color = "#4c4c4c"

        p.line([-self.window, self.window], 
               [1, 1], 
               line_width=self.horizontal.line_width, 
               color=self.horizontal.color, 
               # legend_label=""
               line_dash=self.horizontal.line_dash,
               )
        
        # Add a vertical dotted line at x=0.0
        p.line([0, 0], 
               [ymin, ymax], 
               line_width=self.vertical.line_width, 
               color=self.vertical.color, 
               # legend_label="Vertical Line", 
               line_dash=self.vertical.line_dash,
               )

        title_len = len(self.title)
        width_padding: int = 247
        default_title: int = 11
        # adj_width_padding = width_padding + 6 * abs(title_len - default_title)

        if use_title:
            title_div = Div(
                    text=f"""<div style="background-color: lightgray; width: 100%; height: {self.title_box_height}px; padding: 10px {adj_width_padding}px; 
                      font-size: {self.title_font_size}pt; 
                      font-weight: {self.title_font_weight}; 
                      margin-bottom: 0px; 
                      margin-left: 60px; 
                      text-align: center; 
                      border: {self.title_border_width}px solid black; 
                      display: flex; 
                      justify-content: center;"> """ + f"""</div>""",
                    # width=200,
                    # height=60
            )

        # layout = gridplot([[title_div], [p]], sizing_mode="inherit")
        # layout = column([title_div, p], sizing_mode="inherit")
        if xlabel is None:
            p.xaxis.axis_label = self.xlabel.label 
        else:
            p.xaxis.axis_label = xlabel
        if ylabel is None:
            p.yaxis.axis_label = self.ylabel.label 
        else:
            p.yaxis.axis_label = ylabel
        p.xaxis.axis_label_text_font_style = "normal" 
        p.yaxis.axis_label_text_font_style = "normal" 
        p.xaxis.axis_label_text_font_size = f"{self.xlabel.size}pt"  
        p.yaxis.axis_label_text_font_size = f"{self.ylabel.size}pt" 
        p.xaxis.major_label_text_font_size = f"{self.xlabel.tick_size}pt"
        p.yaxis.major_label_text_font_size = f"{self.ylabel.tick_size}pt"
        self.p = p
        return p
   
def plot_phylum_clustermap(phylum_df):
    pass

def to_png(p, filename: str) -> None:
    export_png(p, filename=filename)
    return

def to_html(p, filename: str) -> None:
    output_file(filename)
    save(p)

if __name__ == "__main__":

    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(description="""Plotting TSS/TES Density data within a given window.""")
    parser.add_argument("--density_data", type=str)
    parser.add_argument("--filename", type=str, default="enrichment_domains.png")
    parser.add_argument("--level", type=str, default="domain")
    parser.add_argument("--mode", type=str, default="density", choices=["density", "template"])

    args = parser.parse_args()
    filename = args.filename
    mode = args.mode
    density_data = args.density_data
    density = pd.read_csv(density_data, index_col=0).drop(columns=['biotype', 'order']).iloc[:3, :]
    level = args.level

    avg = density.iloc[0]
    ci_lower = density.iloc[1]
    ci_upper = density.iloc[2]

    # output_notebook()
    if level == "domain":
        domains = ["Eukaryota", "Bacteria", "Viruses", "Archaea"]
        domain_figures = {}
        for domain in domains:
            plotter = DensityPlotter(title=domain)
            p = plotter.plot(density, ci_lower=ci_lower, ci_upper=ci_upper, mode=mode)
            domain_figures[domain] = p
        figure = gridplot([domain_figures["Eukaryota"], domain_figures["Archaea"]],
                          [domain_figures["Bacteria"], domain_figures["Viruses"]])
    else:
        plotter = DensityPlotter(title=level)
        p = plotter.plot(density, ci_lower=ci_lower, ci_upper=ci_upper, mode=mode)
        figure = p

    show(figure)
    to_png(figure, filename)
