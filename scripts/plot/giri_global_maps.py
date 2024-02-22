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
    """Get all flooded sites
    """
    flood_df = pd.read_csv(os.path.join(output_path,
                "global_exposure_stats_giri",
                "supply_chain_sites_all_climate_rcp_rp_flood_depth.csv"))

    all_flooded = list(set(flood_df.site_id.values.tolist()))
    df = sites[["site_id","site_type","geometry"]]
    df["flooded"] = np.where(df["site_id"].isin(all_flooded),1,0)
    flood_classes = ["not"] + sites_classes
    flood_classes_label = ["not"] + sites_classes_label
    flood_colors = [noflood_color] + class_colors
    
    ax_proj = get_projection(epsg=4326)
    fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(12,6),
                    dpi=500)
    # ax_plots = ax_plots.flatten()
    ax = plot_global_basemap(ax_plots,
                        include_labels=True,
                        label_countries=site_isos)
    z = 5
    legend_handles = []
    for idx,(cl,cc,cb) in enumerate(zip(flood_classes,flood_colors,flood_classes_label)):
        if cl == 'not':
            plot_assets = df[df["flooded"] == 0]
            z_order = z
            color = cc
            fontsize = 4.0
            marker = 'o'
        else:
            plot_assets = df[(df["flooded"] == 1) & (df["site_type"] == cl)]
            z_order = z
            color = cc
            fontsize = 6.0
            marker = 's'
        if len(plot_assets.index) > 0:
            ax = plot_point_assets(ax,plot_assets,
                                    color,
                                    fontsize,
                                    marker,
                                    z_order,cl)
            legend_handles.append(plt.plot([],[],
                                marker=marker, 
                                ms=10.0, 
                                ls="",
                                color=color,
                                label=f"{cb.upper()} FLOODED")[0])
        z += 1
    ax.legend(handles=legend_handles,
                fontsize=10,
                title="$\\bf{Current\ and\ Future\ Site\ Exposures}$",
                loc='lower left')
    plt.tight_layout()
    save_fig(os.path.join(figures,"global_sites_flooded"))
    plt.close()
    
    # """Expected flood depths
    # """
    plot_types = {
                        'climate_scenarios':['baseline','SSP1 RCP2.6','SSP5 RCP8.5'],
                        'climate_scenarios_label': ['Baseline - 2023','SSP1 RCP2.6 - 2080','SSP5 RCP8.5 - 2080'], 
                        'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        'scenario_marker':['s','<','^','v'],     
                    }
    flood_df = pd.read_csv(os.path.join(output_path,
                "global_exposure_stats_giri",
                "giri_supply_chain_sites_all_expected_flood_depths.csv"))

    figure_texts = ['a.','b.','c.']
    flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
    flood_column = "expected_flood_depth_undefended"
    for pt_typ in ["no classes","site_types"]:
        ax_proj = get_projection(epsg=4326)
        fig, ax_plots = plt.subplots(3,1,
                subplot_kw={'projection': ax_proj},
                figsize=(12,16),
                dpi=500)
        ax_plots = ax_plots.flatten()
        j = 0
        for idx,(sc,sb) in enumerate(zip(plot_types['climate_scenarios'],
                        plot_types['climate_scenarios_label'])):
            df = flood_df[flood_df["rcp"].isin([sc])]
            values_range = flood_df[flood_column].values.tolist()
            df = pd.merge(sites[[asset_id,"site_type","geometry"]],
                            df[[asset_id,flood_column]],
                            how="left",on=[asset_id]).fillna(0)
            legend_handles = []
            ax = plot_global_basemap(ax_plots[j],
                            include_labels=
                            True,label_countries=site_isos,
                            scalebar_location=(0.60,0.05),
                            arrow_location=(0.54,0.08))
            
            if pt_typ == "no classes":
                ax = point_map_plotting_colors_width(ax,df,flood_column,values_range,
                            point_colors=flood_colors,
                            point_zorder=[6,7,8,9,10],
                            legend_label="Expected Annual flood depth (m)",
                            no_value_label="Site not flooded",
                            no_value_color="#969696",
                            point_steps=len(flood_colors) + 1,
                            width_step=9.0,
                            interpolation="linear",
                            legend_size=12.0,
                            legend_weight=2.0
                            )
            else:
                ax = point_map_plotting_colors_width(ax,df,
                            flood_column,
                            values_range,
                            point_classify_column="site_type",
                            point_categories=sites_classes,
                            point_colors=class_colors,
                            point_labels=[s.upper() for s in sites_classes_label],
                            point_zorder=[6,7,8,9],
                            width_step = 9.0,
                            interpolation = 'linear',
                            legend_label="Expected Annual flood depth (m)",
                            legend_size=12,
                            legend_weight=2.0,
                            no_value_label="Site not flooded",
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
            f"global_expected_flood_depths_{pt_typ}.png"))
        plt.close()

    """Expected flood depths China, Japan, Korea
    """
    plot_types = {
                        'climate_scenarios':['baseline','SSP1 RCP2.6','SSP5 RCP8.5'],
                        'climate_scenarios_label': ['Baseline - 2023','SSP1 RCP2.6 - 2080','SSP5 RCP8.5 - 2080'], 
                        'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        'scenario_marker':['s','<','^','v'],     
                    }
    flood_df = pd.read_csv(os.path.join(output_path,
                "global_exposure_stats_giri",
                "giri_supply_chain_sites_all_expected_flood_depths.csv"))
    chn_jpn_kor_sites = sites[sites["iso_code"].isin(["CHN","JPN","KOR","PRK"])]

    figure_texts = ['a.','b.','c.']
    flood_colors = ['#c994c7','#4eb3d3','#0570b0','#6a51a3','#4d004b']
    flood_column = "expected_flood_depth_undefended"
    for pt_typ in ["no classes","site_types"]:
        ax_proj = get_projection(epsg=4326)
        fig, ax_plots = plt.subplots(3,1,
                subplot_kw={'projection': ax_proj},
                figsize=(12,16),
                dpi=500)
        ax_plots = ax_plots.flatten()
        j = 0
        for idx,(sc,sb) in enumerate(zip(plot_types['climate_scenarios'],
                        plot_types['climate_scenarios_label'])):
            df = flood_df[flood_df["rcp"].isin([sc])]
            values_range = flood_df[flood_column].values.tolist()
            df = pd.merge(chn_jpn_kor_sites[[asset_id,"site_type","geometry"]],
                            df[[asset_id,flood_column]],
                            how="left",on=[asset_id]).fillna(0)
            legend_handles = []
            ax = plot_global_basemap(ax_plots[j],
                                include_countries=["CHN","JPN","KOR","PRK"],
                                include_labels=True,label_countries=site_isos,
                                scalebar_location=(0.60,0.05),
                                arrow_location=(0.54,0.08),label_size=15)
            
            if pt_typ == "no classes":
                ax = point_map_plotting_colors_width(ax,df,flood_column,values_range,
                            point_colors=flood_colors,
                            point_zorder=[6,7,8,9,10],
                            legend_label="Expected Annual flood depth (m)",
                            no_value_label="Site not flooded",
                            no_value_color="#969696",
                            point_steps=len(flood_colors) + 1,
                            width_step=9.0,
                            interpolation="linear",
                            legend_size=12.0,
                            legend_weight=2.0
                            )
            else:
                ax = point_map_plotting_colors_width(ax,df,
                            flood_column,
                            values_range,
                            point_classify_column="site_type",
                            point_categories=sites_classes,
                            point_colors=class_colors,
                            point_labels=[s.upper() for s in sites_classes_label],
                            point_zorder=[6,7,8,9],
                            width_step = 9.0,
                            interpolation = 'linear',
                            legend_label="Expected Annual flood depth (m)",
                            legend_size=12,
                            legend_weight=2.0,
                            no_value_label="Site not flooded",
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
            f"chn_jpn_kor_expected_flood_depths_{pt_typ}.png"))
        plt.close()


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
