"""GIRI flood map plots
"""
import os
import sys
from collections import OrderedDict
import pandas as pd
import geopandas as gpd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from map_plotting_utils import (load_config,get_projection,
                            plot_basemap,plot_global_basemap, plot_point_assets, 
                            point_map_plotting_color_width,
                            point_map_plotting_colors_width, save_fig)
from tqdm import tqdm
tqdm.pandas()


def main(config):
    data_path = config['paths']['data']
    output_path = config['paths']['output']
    figure_path = config['paths']['figures']

    figures = os.path.join(figure_path,"global_flood_maps")
    if os.path.exists(figures) == False:
        os.mkdir(figures)

    
    """Plot of sample global flood map
    """
    flood_file = os.path.join(processed_data_path,
                            "incoming_data",
                            "unepgrid-2023-giri",
                            "global_pc_h100glob.tif")
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    # ax_plots = ax_plots.flatten()
    ax = plot_global_basemap(ax_plots)
    im = plot_raster(ax, flood_file,cmap="terrain")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
    # cbar.set_clim(vmin=0,vmax=max_val)

    cbar.outline.set_color("none")
    cbar.ax.yaxis.set_tick_params(color='black')
    cbar.ax.set_xlabel('Flood depths (m)',fontsize=7,color='black')

    plt.title("Baseline - 1 in 100 year river flooding", fontsize = 10)
    save_fig(os.path.join(figures,"global_river_flood_map_100_year_return_period.png"))


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
