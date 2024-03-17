"""GIRI exposure analysis results
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

    figures = os.path.join(figure_path,"global_supply_chain_sites_giri")
    sites = gpd.read_file(os.path.join(
                            data_path,
                            'supply_chain_study',
                            'supply_chain_sites_all_qingnan_edited.gpkg'),
                            layer="nodes")
    site_isos = list(set(sites.iso_code.values.tolist()))
    sites_classes = ["buyer",
                        "3 major producers in China",
                        "8 other producers in China",
                        "alternative producers owned by 3 non-Chinese firms"
                    ]
    sites_classes_label = ["buyer",
                        "3 major producers in China",
                        "8 other producers in China",
                        "producers owned by 3 non-Chinese firms"
                    ]
    class_colors = ["#3182bd",
                        "#de2d26",
                        "#67000d",
                        '#31a354']
    noflood_color = '#969696'
    asset_id = "site_id"
    
    # """Expected flood depths
    # """
    flood_df = pd.read_csv(os.path.join(output_path,
                "global_exposure_stats_giri",
                "giri_supply_chain_sites_all_expected_flood_depths.csv"))
    flood_column = "expected_flood_depth_undefended"
    f_df = pd.DataFrame(list(set(flood_df[asset_id].values.tolist())),columns=[asset_id])
    cl_sc = ['baseline','SSP1 RCP2.6','SSP5 RCP8.5']
    cl_sc_col = ['baseline','rcp26','rcp85']
    for idx,(cs,cn) in enumerate(zip(cl_sc,cl_sc_col)):
        fdf = flood_df[flood_df["rcp"] == cs]
        fdf.rename(columns={flood_column:cn},inplace=True)
        f_df = pd.merge(f_df,fdf[[asset_id,cn]],how="left",on=[asset_id]).fillna(0)

    f_df["rcp26_diff"] = f_df["rcp26"] - f_df["baseline"]
    f_df["rcp85_diff"] = f_df["rcp85"] - f_df["baseline"]
    f_df.to_csv(os.path.join(output_path,
                "global_exposure_stats_giri",
                "giri_supply_chain_sites_all_expected_flood_depth_changes.csv"),index=False)


    plot_types = {
                        'climate_scenarios':['baseline','rcp26_diff','rcp85_diff'],
                        'climate_scenarios_label': ['Baseline - 2023','SSP1 RCP2.6 - 2080','SSP5 RCP8.5 - 2080'], 
                        'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        'scenario_marker':['s','<','^','v'],     
                    }
    figure_texts = ['a.','b.','c.']
    flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(3,1,
            subplot_kw={'projection': ax_proj},
            figsize=(12,16),
            dpi=500)
    ax_plots = ax_plots.flatten()
    j = 0
    for idx,(sc,sb) in enumerate(zip(plot_types['climate_scenarios'],
                    plot_types['climate_scenarios_label'])):
        df = pd.merge(sites[[asset_id,"site_type","geometry"]],
                        f_df[[asset_id,sc]],
                        how="left",on=[asset_id]).fillna(0)
        legend_handles = []
        ax = plot_global_basemap(ax_plots[j],
                        include_labels=
                        True,label_countries=site_isos,
                        scalebar_location=(0.60,0.05),
                        arrow_location=(0.54,0.08))
        values_range = df[sc].values.tolist()
        if sc == "baseline":
            label = "Expected Annual flood depth (m)"
            flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
        else:
            label = "Expected Annual flood depth change (m)"
            flood_colors = ['#1a9850','#a6d96a','#fee08b','#fdae61','#d73027']

        ax = point_map_plotting_colors_width(ax,df,sc,values_range,
                    point_colors=flood_colors,
                    point_zorder=[6,7,8,9,10],
                    legend_label=label,
                    no_value_label="Site not flooded",
                    no_value_color="#969696",
                    point_steps=len(flood_colors) + 1,
                    width_step=9.0,
                    interpolation="linear",
                    legend_size=12.0,
                    legend_weight=2.0
                    )
        
        text_label = f"{figure_texts[j]} {sb}"
        ax.text(
        0.05,
        0.95,
        text_label,
        horizontalalignment='left',
        transform=ax.transAxes,
        size=18,
        weight='bold')     

        j+=1            

    plt.tight_layout()
    save_fig(os.path.join(figures,
        f"global_expected_flood_depths_changes.png"))
    plt.close()

    """Expected flood depths China, Japan, Korea
    """
    plot_types = {
                        'climate_scenarios':['baseline','rcp26_diff','rcp85_diff'],
                        'climate_scenarios_label': ['Baseline - 2023','SSP1 RCP2.6 - 2080','SSP5 RCP8.5 - 2080'], 
                        'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        'scenario_marker':['s','<','^','v'],     
                    }
    chn_jpn_kor_sites = sites[sites["iso_code"].isin(["CHN","JPN","KOR","PRK"])]

    figure_texts = ['a.','b.','c.']
    flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(3,1,
            subplot_kw={'projection': ax_proj},
            figsize=(12,16),
            dpi=500)
    ax_plots = ax_plots.flatten()
    j = 0
    for idx,(sc,sb) in enumerate(zip(plot_types['climate_scenarios'],
                    plot_types['climate_scenarios_label'])):
        df = pd.merge(chn_jpn_kor_sites[[asset_id,"site_type","geometry"]],
                        f_df[[asset_id,sc]],
                        how="left",on=[asset_id]).fillna(0)
        values_range = df[sc].values.tolist()
        legend_handles = []
        ax = plot_global_basemap(ax_plots[j],
                            include_countries=["CHN","JPN","KOR","PRK"],
                            include_labels=True,label_countries=site_isos,
                            scalebar_location=(0.60,0.05),
                            arrow_location=(0.54,0.08),label_size=15)
        
        if sc == "baseline":
            label = "Expected Annual flood depth (m)"
            flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
        else:
            label = "Expected Annual flood depth change (m)"
            flood_colors = ['#1a9850','#a6d96a','#fee08b','#fdae61','#d73027']
        ax = point_map_plotting_colors_width(ax,df,sc,values_range,
                    point_colors=flood_colors,
                    point_zorder=[6,7,8,9,10],
                    legend_label=label,
                    no_value_label="Site not flooded",
                    no_value_color="#969696",
                    point_steps=len(flood_colors) + 1,
                    width_step=9.0,
                    interpolation="linear",
                    legend_size=12.0,
                    legend_weight=2.0
                    )
        text_label = f"{figure_texts[j]} {sb}"
        ax.text(
        0.05,
        0.95,
        text_label,
        horizontalalignment='left',
        transform=ax.transAxes,
        size=18,
        weight='bold')     

        j+=1            

    plt.tight_layout()
    save_fig(os.path.join(figures,
        f"chn_jpn_kor_expected_flood_depths_changes.png"))
    plt.close()


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
