"""Road network risks and adaptation maps
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
                            plot_basemap, plot_point_assets, 
                            plot_line_assets, save_fig)
from tqdm import tqdm
tqdm.pandas()

def main(config):
    data_path = config['paths']['data']
    output_path = config['paths']['output']
    figure_path = config['paths']['figures']

    threshold = 0.0
    threshold_label = 'buyer'
    exposure_results_path = os.path.join(output_path,
                f'global_exposure_stats') # Where we have all the risk results

    figures = os.path.join(figure_path,"global_supply_chain_sites")
    if os.path.exists(figures) == False:
        os.mkdir(figures)

    asset_id = "site_id"
    boundary_id = "province"
    # company_id = "Company_ENG"
    company_id = "company_name_eng_merged"

    plot_types = [
                        {
                        'type':'all_exposures',
                        'percentage':False,
                        'groupby':[
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'years':[2030,2050,2080],
                        'climate_scenarios':['rcp4p5','rcp8p5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_rcp_rp_epoch.xlsx',
                        'plot_name':'numbers_by_rcp_rp_epoch.png'
                        },
                        {
                        'type':'all_exposures',
                        'percentage':True,
                        'groupby':[
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'years':[2030,2050,2080],
                        'climate_scenarios':['rcp4p5','rcp8p5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_rcp_rp_epoch.xlsx',
                        'plot_name':'percentage_by_rcp_rp_epoch.png'
                        },
                        {
                        'type':'company_type',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'years':[2030,2050,2080],
                        'climate_scenarios':['rcp4p5','rcp8p5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_site_type_rcp_rp_epoch.xlsx',
                        'plot_name':'numbers_by_rcp_rp_epoch.png'
                        },
                        {
                        'type':'company_type',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'years':[2030,2050,2080],
                        'climate_scenarios':['rcp4p5','rcp8p5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_site_type_rcp_rp_epoch.xlsx',
                        'plot_name':'percentage_by_rcp_rp_epoch.png'
                        },
                        {
                        'type':'company_country',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'country',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'site_groupby':[
                                    'site_type',
                                    'country'
                                ],
                        'x_groupby':'country',
                        'sites_classes':["buyer",
                                    "3 major producers in China",
                                    "8 other producers in China",
                                    "alternative producers owned by 3 non-Chinese firms"
                                ],
                        'sites_classes_label':["buyer",
                                            "3 major producers in China",
                                            "8 other producers in China",
                                            "producers owned by 3 non-Chinese firms"
                                        ],
                        'class_colors':["#3182bd",
                                            "#de2d26",
                                            "#67000d",
                                            '#31a354'],    
                        'file_name':'exposures_numbers_by_company_country_rcp_rp_epoch.xlsx',
                        'plot_name':'max_numbers_across_hazards.png'
                        },
                        {
                        'type':'company_country',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'country',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'site_groupby':[
                                    'site_type',
                                    'country'
                                ],
                        'x_groupby':'country',
                        'sites_classes':["buyer",
                                    "3 major producers in China",
                                    "8 other producers in China",
                                    "alternative producers owned by 3 non-Chinese firms"
                                ],
                        'sites_classes_label':["buyer",
                                            "3 major producers in China",
                                            "8 other producers in China",
                                            "producers owned by 3 non-Chinese firms"
                                        ],
                        'class_colors':["#3182bd",
                                            "#de2d26",
                                            "#67000d",
                                            '#31a354'], 
                        'file_name':'exposures_numbers_by_company_country_rcp_rp_epoch.xlsx',
                        'plot_name':'max_percentages_across_hazards.png'
                        },
                        {
                        'type':'company_name',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'company_name_eng_merged',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'site_groupby':[
                                    'site_type',
                                    'company_name_eng_merged'
                                ],
                        'x_groupby':'company_name_eng_merged',        
                        'sites_classes':["buyer",
                                    "3 major producers in China",
                                    "8 other producers in China",
                                    "alternative producers owned by 3 non-Chinese firms"
                                ],
                        'sites_classes_label':["buyer",
                                            "3 major producers in China",
                                            "8 other producers in China",
                                            "producers owned by 3 non-Chinese firms"
                                        ],
                        'class_colors':["#3182bd",
                                            "#de2d26",
                                            "#67000d",
                                            '#31a354'], 
                        'file_name':'exposures_numbers_by_company_rcp_rp_epoch.xlsx',
                        'plot_name':'company_name_max_numbers_across_hazards.png'
                        },
                        {
                        'type':'company_name',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'company_name_eng_merged',
                                    'rcp',
                                    'epoch',
                                    'rp'
                                ],
                        'site_groupby':[
                                    'site_type',
                                    'company_name_eng_merged'
                                ],
                        'x_groupby':'company_name_eng_merged',
                        'sites_classes':["buyer",
                                    "3 major producers in China",
                                    "8 other producers in China",
                                    "alternative producers owned by 3 non-Chinese firms"
                                ],
                        'sites_classes_label':["buyer",
                                            "3 major producers in China",
                                            "8 other producers in China",
                                            "producers owned by 3 non-Chinese firms"
                                        ],
                        'class_colors':["#3182bd",
                                            "#de2d26",
                                            "#67000d",
                                            '#31a354'], 
                        'file_name':'exposures_numbers_by_company_rcp_rp_epoch.xlsx',
                        'plot_name':'company_name_max_percentages_across_hazards.png'
                        },
                        # {
                        # 'type':'exposures_across_all_hazards',
                        # 'groupby':[],
                        # 'years':[],
                        # 'climate_scenarios':[], 
                        # 'scenario_color':[], 
                        # 'scenario_marker':[],     
                        # # 'file_name':'supply_chain_sites_lithium_units_splits__aqueduct_river__nodes.csv',
                        # 'file_name':'supply_chain_sites_all_splits__aqueduct_river__nodes.csv',
                        # 'plot_name':'exposures_across_all_hazards.png'
                        # },
                        # {
                        # 'type':'exposures_over_time',
                        # 'groupby':[],
                        # 'years':[1980,2030,2050,2080],
                        # 'climate_scenarios':['rcp4p5','rcp8p5'], 
                        # 'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        # 'scenario_marker':['s','<','^','v'],     
                        # # 'file_name':'supply_chain_sites_lithium_units_splits__aqueduct_river__nodes.csv',
                        # 'file_name':'supply_chain_sites_all_splits__aqueduct_river__nodes.csv',
                        # 'plot_name':'exposures_over_time.png'
                        # },
                        # {
                        # 'type':'exposures_climate_scenarios_time_epochs',
                        # 'groupby':[],
                        # 'years':[1980,2030,2050,2080],
                        # 'climate_scenarios':['rcp4p5','rcp8p5'], 
                        # 'scenario_color':["#9ecae1","#2171b5","#08306b","#54278f"], 
                        # 'scenario_marker':['s','<','^','v'],     
                        # # 'file_name':'supply_chain_sites_lithium_units_splits__aqueduct_river__nodes.csv',
                        # 'file_name':'supply_chain_sites_all_splits__aqueduct_river__nodes.csv',
                        # 'plot_name':'exposures_climate_scenarios_time_epochs.png'
                        # },
                    ]
    china_codes = ["MAC","HKG","TWN","CHN"]
    # sites = gpd.read_file(os.path.join(data_path,
    #                         "supply_chain_study",
    #                         "supply_chain_sites_lithium_units.gpkg"),
    #                     layer="nodes")
    sites = gpd.read_file(os.path.join(
                            data_path,
                            'supply_chain_study',
                            'supply_chain_sites_all_qingnan_edited.gpkg'),
                            layer="nodes")
    # total_china_sites = len(sites[sites["iso_code"].isin(china_codes)])
    total_china_sites = len(sites.index)
    baseyear = 2023
    flood_color = '#3182bd'
    noflood_color = '#969696' 
    protection_standards_colors = [(0,'#feb24c'),
                                (10,'#fd8d3c'),
                                (20,'#fc4e2a'),
                                (25,'#e31a1c'),
                                (50,'#bd0026'),
                                (100,'#800026'),
                                ('Not flooded',noflood_color)]
    boundary_gpd = gpd.read_file(os.path.join(data_path,
                                'admin_boundaries',
                                'China_regions.gpkg'),encoding="utf-8")
    bounds = boundary_gpd.geometry.total_bounds # this gives your boundaries of the map as (xmin,ymin,xmax,ymax)
    del boundary_gpd

    hazard_data_details = pd.read_csv(os.path.join(data_path,
                                    "aqueduct_river.csv"),encoding="latin1")
    hazard_keys = hazard_data_details["key"].values.tolist()
    """Get all the exposure plots
    """
    for st in range(len(plot_types)):
        if plot_types[st]["type"] in ("all_exposures","company_type"):
            quantiles_list = ['mean','min','max']
            exposures_df = pd.read_excel(os.path.join(exposure_results_path,
                                            f"{plot_types[st]['file_name']}"),
                                sheet_name="exposures")[
                                plot_types[st]['groupby'] + [f"flood_count_{g}" for g in quantiles_list]]
            if plot_types[st]["type"] == "all_exposures":
                pt_loop = [plot_types[st]["type"]]
            else:
                pt_loop = list(set(exposures_df["site_type"].values.tolist()))

            for i in pt_loop:
                if i != plot_types[st]["type"]:
                    exposures = exposures_df[exposures_df["site_type"] == i]
                    site_numbers = len(sites[sites['site_type'] == i].index)
                else:
                    exposures = exposures_df.copy()
                    site_numbers = len(sites.index)
                if plot_types[st]["percentage"] is False:
                    factor = 1
                    ylabel = "Flooded number of sites"
                else:
                    factor = 100.0/site_numbers
                    ylabel  = "Flooded percentage of sites (%)"

                figure_texts = ['a.','b.','c.']
                fig, ax_plots = plt.subplots(1,3,
                        figsize=(30,10),
                        dpi=500)
                ax_plots = ax_plots.flatten()
                j = 0
                for year in plot_types[st]['years']:
                    # exposures = exposures[exposures["rp"] > 2.0]
                    y_min = factor*min(exposures[[f"flood_count_{g}" for g in quantiles_list]].min())
                    y_max = factor*max(exposures[[f"flood_count_{g}" for g in quantiles_list]].max())
                    rps = list(set(exposures['rp'].values)) 
                    ax = ax_plots[j]
                    ax.plot(exposures[exposures['epoch'] == baseyear]['rp'],
                            factor*exposures[exposures['epoch'] == baseyear][f"flood_count_{quantiles_list[0]}"],
                            'o-',color='#fd8d3c',markersize=10,linewidth=2.0,
                            label='Baseline')
                    for c in range(len(plot_types[st]['climate_scenarios'])):
                        sc = plot_types[st]['climate_scenarios'][c]
                        cl = plot_types[st]['scenario_color'][c]
                        m = plot_types[st]['scenario_marker'][c]
                        exp = exposures[(exposures['epoch'] == year) & (exposures['rcp'] == sc)]
                        ax.plot(exp['rp'],
                                factor*exp[f"flood_count_{quantiles_list[0]}"],
                                m,color=cl,markersize=10,linewidth=2.0,
                                label=f"RCP {sc.upper()} mean")
                        # ax.fill_between(exp['rp'],sector['length_unit']*exp[f"{sector['exposure_column']}_{quantiles_list[0]}"],
                        #     sector['length_unit']*exp[f"{sector['exposure_column']}_{quantiles_list[1]}"],alpha=0.3,facecolor=cl)
                        ax.fill_between(exp['rp'],factor*exp[f"flood_count_{quantiles_list[1]}"],
                            factor*exp[f"flood_count_{quantiles_list[2]}"],alpha=0.3,facecolor=cl,
                            label=f"RCP {sc.upper()} min-max")


                    
                    ax.set_xlabel('Return period (years)',fontsize=16,fontweight='bold')
                    ax.set_ylabel(ylabel,fontsize=16,fontweight='bold')
                    ax.set_xscale('log')
                    ax.set_ylim([y_min,y_max])
                    ax.tick_params(axis='both', labelsize=14)
                    ax.set_xticks([t for t in rps])
                    ax.set_xticklabels([str(t) for t in rps])
                    ax.grid(True)
                    # ax.set_xticks([t for t in list(set(exposures[exposures['year'] == baseyear]['return_period'].values))], 
                    #             [str(t) for t in list(set(exposures[exposures['year'] == baseyear]['return_period'].values))])
                    ax.text(
                        0.05,
                        0.95,
                        figure_texts[j],
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=18,
                        weight='bold')
                    ax.text(
                        0.05,
                        0.80,
                        year,
                        horizontalalignment='left',
                        transform=ax.transAxes,
                        size=18,
                        weight='bold')
                    # if j <= 3:
                    #     # ax.text(
                    #     #     0.35,
                    #     #     0.95,
                    #     #     sector['sector_label'],
                    #     #     horizontalalignment='left',
                    #     #     transform=ax.transAxes,
                    #     #     size=18,
                    #     #     weight='bold')
                    #     ax.set_title(
                    #             sector['sector_label'],
                    #             fontdict = {'fontsize':18,
                    #             'fontweight':'bold'})

                    j+=1            

                ax_plots[2].legend(
                            loc='upper left', 
                            bbox_to_anchor=(1.05,0.8),
                            prop={'size':18,'weight':'bold'})
                plt.tight_layout()
                save_fig(os.path.join(figures,f"{i}_{plot_types[st]['plot_name']}"))
                plt.close()
        
        elif plot_types[st]['type'] in ("company_country","company_name"):
            exposures = pd.read_excel(os.path.join(exposure_results_path,
                                    f"{plot_types[st]['file_name']}"))
            exposures = exposures.groupby(plot_types[st]['groupby'])["flood_count_max"].sum().reset_index()
            site_numbers = sites.groupby(plot_types[st]['site_groupby']).size().reset_index(name='counts')
            exposures = pd.merge(exposures,site_numbers,how="left",on=plot_types[st]['site_groupby'])
            exposures["percentage_flooded"] = 100.0*exposures["flood_count_max"]/exposures["counts"]
            if "country" in exposures.columns.values.tolist():
                exposures["country"] = exposures["country"].str.replace("of America","")
            if plot_types[st]["percentage"] is False:
                ylabel = "Flooded number of sites"
                exposure_column = "flood_count_max"
            else:
                ylabel  = "Flooded percentage of each site type (%)"
                exposure_column = "percentage_flooded"

            exposures = exposures.sort_values(by=exposure_column,ascending=False)
            exposures = exposures.drop_duplicates(subset=plot_types[st]['site_groupby'],keep="first") 
            country_totals = exposures.groupby(plot_types[st]['x_groupby'])[exposure_column].sum().reset_index()
            country_totals = country_totals.sort_values(by=exposure_column,ascending=False)
            ymax = 1.05*max(country_totals[exposure_column])
            all_countries = country_totals[plot_types[st]['x_groupby']].values.tolist()
            fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=500)
            bottom  = np.zeros(len(all_countries))
            for idx,(site_type,site_type_color,site_label) in enumerate(list(zip(plot_types[st]['sites_classes'],
                                                            plot_types[st]['class_colors'],
                                                            plot_types[st]['sites_classes_label']))): 
                df = exposures[exposures.site_type == site_type]
                if len(df.index) > 0:
                    vals = []
                    for hz in all_countries:
                        df1 = df[df[plot_types[st]['x_groupby']] == hz]
                        if len(df1.index) > 0:
                            vals.append(df1[exposure_column].values[0])
                        else:
                            vals.append(0)

                    vals = np.array(vals)
                    ax.bar(all_countries,vals,
                                width=0.5,
                                color=site_type_color,
                                alpha=0.5, 
                                label=f"{site_label.upper()}",
                                bottom =bottom)
                    bottom += vals

                for c in ax.containers:
                    # Optional: if the segment is small or 0, customize the labels
                    # labels = [round (v.get_height(),0) if v.get_height() > 5.0 else '' for v in c]
                    labels = [int(v.get_height()) if v.get_height() > 0.0 else '' for v in c]
                    
                    # remove the labels parameter if it's not needed for customized labels
                    ax.bar_label(c, labels=labels, label_type='center')

            # ax.legend(prop={'size':12,'weight':'bold'},bbox_to_anchor=(1.5, 0.95))
            ax.legend(prop={'size':12,'weight':'bold'})
            ax.set_ylim([0,1.2*ymax])
            ax.set_ylabel(ylabel,fontweight='bold',fontsize=15)
            # xticks = list(np.arange(1,len(all_hazards),1, dtype=int))
            if plot_types[st]['x_groupby'] == "company_name_eng_merged":
                all_countries = [a.replace(" INDUSTRY CO.,LTD. (Tyeeli)","") for a in all_countries]
                ac = [str(ac[:3]).upper() for ac in all_countries]
                dcsummary = [f"{a} - {b}" for (a,b) in zip(ac,all_countries)]
                dcsummary = np.reshape(dcsummary, (7, 3))
                the_table = plt.table(cellText=dcsummary,
                        colWidths = [0.4]*3,
                        cellLoc = 'center', rowLoc = 'center',
                        loc='bottom', bbox=[-0.1, -0.4, 1.2, 0.3])
                the_table.auto_set_font_size(False)
                the_table.set_fontsize(9)
                all_countries = ac
            
            ax.set_xticklabels(all_countries,
                            fontsize=9, rotation=45,
                            fontweight="bold")
            # ax.set_axisbelow(True)
            # ax.grid(True)
            plt.tight_layout()
            save_fig(os.path.join(figures,plot_types[st]['plot_name']))
            plt.close()
        elif plot_types[st]['type'] == "exposures_across_all_hazards":
            exposures = pd.read_csv(os.path.join(exposure_results_path,
                                    f"{plot_types[st]['file_name']}"))
            exposures["all_exposures"] = exposures[hazard_keys].sum(axis=1) 
            exposures = pd.merge(exposures[[asset_id,"all_exposures"]],sites,how="left",on=[asset_id])
            exposures = gpd.GeoDataFrame(exposures[exposures["iso_code"].isin(china_codes)],
                                    geometry="geometry",crs="EPSG:4326")
            ax_proj = get_projection(extent = (bounds[0]+5,bounds[2]-10,bounds[1],bounds[3]))
            fig, ax_plots = plt.subplots(1,1,
                    subplot_kw={'projection': ax_proj},
                    figsize=(10,8),
                    dpi=500)
            # ax_plots = ax_plots.flatten()
            ax = plot_basemap(ax_plots,include_labels=True)
            legend_handles = []
            for psc in ["Flooded","Not flooded"]:
                if psc == 'Not flooded':
                    plot_assets = exposures[exposures["all_exposures"] == 0]
                    z_order = 7
                    color = noflood_color
                    marker = 'o'
                else:
                    plot_assets = exposures[exposures["all_exposures"] > 0]
                    z_order = 8
                    color = flood_color
                    marker = 's'
                if len(plot_assets.index) > 0:
                    ax = plot_point_assets(ax,plot_assets,
                                            color,
                                            20.0,
                                            marker,
                                            z_order)
                    legend_handles.append(plt.plot([],[],
                                        marker=marker, 
                                        ms=12.0, 
                                        ls="",
                                        color=color,
                                        label=psc)[0])
            ax.legend(handles=legend_handles,fontsize=10,title="$\\bf{Current\ and\ Future\ Site\ Exposures}$",loc='lower left')
            # ax.text(
            #     0.05,
            #     0.95,
            #     figure_texts[j+4],
            #     horizontalalignment='left',
            #     transform=ax.transAxes,
            #     size=18,
            #     weight='bold')     

            # j+=1            

            plt.tight_layout()
            save_fig(os.path.join(figures,plot_types[st]['plot_name']))
            plt.close()
        elif plot_types[st]['type'] == "exposures_over_time":
            figure_texts = ['a.','b.']
            ax_proj = get_projection(extent = (bounds[0]+5,bounds[2]-10,bounds[1],bounds[3]))
            fig, ax_plots = plt.subplots(1,2,
                    subplot_kw={'projection': ax_proj},
                    figsize=(24,10),
                    dpi=500)
            ax_plots = ax_plots.flatten()
            exposures = pd.read_csv(os.path.join(exposure_results_path,
                                    f"{plot_types[st]['file_name']}"))
            china_ids = sites[sites["iso_code"].isin(china_codes)][asset_id].values.tolist()
            j = 0
            for sc in plot_types[st]['climate_scenarios']:
                legend_handles = []
                ax = plot_basemap(ax_plots[j],include_labels=True)
                hazards = hazard_data_details[hazard_data_details["rcp"].isin(["historical",sc])]
                year_ids = []
                for y in range(len(plot_types[st]['years'])):
                    year = plot_types[st]['years'][y]
                    years_keys = hazards[hazards["epoch"] == year]["key"].values.tolist() 
                    exposures["all_exposures"] = exposures[years_keys].sum(axis=1) 
                    flood_ids = exposures[exposures["all_exposures"]>0][asset_id].values.tolist()
                    year_ids += flood_ids
                    if flood_ids != year_ids:
                        flood_ids = list(set(year_ids)^set(flood_ids))
                    plot_assets = sites[sites[asset_id].isin(flood_ids)]
                    z_order = 8 + j
                    color = plot_types[st]['scenario_color'][y]
                    marker = plot_types[st]['scenario_marker'][y]
                    if year < 2080:
                        label = f"Sites flooded from {year}-2080"
                    else:
                        label = f"Sites flooded from 2080 onwards"
                    if len(plot_assets.index) > 0:
                        ax = plot_point_assets(ax,plot_assets,
                                                color,
                                                30.0,
                                                marker,
                                                z_order)
                        legend_handles.append(plt.plot([],[],
                                            marker=marker, 
                                            ms=12.0, 
                                            ls="",
                                            color=color,
                                            label=label)[0])
                remaining_ids = list(set(china_ids)^set(year_ids))
                if len(remaining_ids) > 0:
                    plot_assets = sites[sites[asset_id].isin(remaining_ids)]
                    ax = plot_point_assets(ax,plot_assets,
                                            noflood_color,
                                            24.0,
                                            'o',
                                            7)
                    legend_handles.append(plt.plot([],[],
                                        marker='o', 
                                        ms=12.0, 
                                        ls="",
                                        color=noflood_color,
                                        label="Not flooded")[0])
                ax.legend(handles=legend_handles,fontsize=10,title="$\\bf{Current\ and\ Future\ Site\ Exposures}$",loc='lower left')
                ax.text(
                    0.05,
                    0.95,
                    f"{figure_texts[j]} RCP {sc}",
                    horizontalalignment='left',
                    transform=ax.transAxes,
                    size=18,
                    weight='bold')     

                j+=1            

            plt.tight_layout()
            save_fig(os.path.join(figures,plot_types[st]['plot_name']))
            plt.close()
        elif plot_types[st]['type'] == "exposures_climate_scenarios_time_epochs":
            figure_texts = ['a.','b.','c.','d.','e.','f.','g.','h.']
            ax_proj = get_projection(extent = (bounds[0]+5,bounds[2]-10,bounds[1],bounds[3]))
            fig, ax_plots = plt.subplots(2,4,
                    subplot_kw={'projection': ax_proj},
                    figsize=(24,10),
                    dpi=500)
            ax_plots = ax_plots.flatten()
            exposures = pd.read_csv(os.path.join(exposure_results_path,
                                    f"{plot_types[st]['file_name']}"))
            # china_ids = sites[sites["iso_code"].isin(china_codes)][asset_id].values.tolist()
            j = 0
            for sc in plot_types[st]['climate_scenarios']:
                hazards = hazard_data_details[hazard_data_details["rcp"].isin(["historical",sc])]
                for y in range(len(plot_types[st]['years'])):
                    legend_handles = []
                    ax = plot_basemap(ax_plots[j],include_labels=True)
                    year = plot_types[st]['years'][y]
                    years_keys = hazards[hazards["epoch"] == year]["key"].values.tolist() 
                    exposures["all_exposures"] = exposures[years_keys].sum(axis=1) 
                    assets = pd.merge(exposures[[asset_id,"all_exposures"]],sites,how="left",on=[asset_id])
                    assets = gpd.GeoDataFrame(assets[assets["iso_code"].isin(china_codes)],
                                                    geometry="geometry",crs="EPSG:4326")
                    for psc in ["Flooded","Not flooded"]:
                        if psc == 'Not flooded':
                            plot_assets = assets[assets["all_exposures"] == 0]
                            z_order = 7
                            color = noflood_color
                            marker = 'o'
                        else:
                            plot_assets = assets[assets["all_exposures"] > 0]
                            z_order = 8
                            color = plot_types[st]['scenario_color'][y]
                            marker = plot_types[st]['scenario_marker'][y]
                        if len(plot_assets.index) > 0:
                            ax = plot_point_assets(ax,plot_assets,
                                                    color,
                                                    20.0,
                                                    marker,
                                                    z_order)
                            legend_handles.append(plt.plot([],[],
                                                marker=marker, 
                                                ms=12.0, 
                                                ls="",
                                                color=color,
                                                label=psc)[0])
                    ax.legend(handles=legend_handles,fontsize=10,title="$\\bf{Site\ Exposures}$",loc='lower left')
                    if year == 1980:
                        text_label = f"{figure_texts[j]} Baseline - {year}"
                    else:
                        text_label = f"{figure_texts[j]} RCP {sc} - {year}"
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
            save_fig(os.path.join(figures,plot_types[st]['plot_name']))
            plt.close()



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
