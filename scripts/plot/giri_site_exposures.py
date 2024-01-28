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
                'global_exposure_stats_giri') # Where we have all the risk results

    figures = os.path.join(figure_path,"global_supply_chain_sites_giri")
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
                                    'rp'
                                ],
                        'climate_scenarios':['SSP1 RCP2.6','SSP5 RCP8.5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_rcp_rp.xlsx',
                        'plot_name':'numbers_by_rcp_rp.png'
                        },
                        {
                        'type':'all_exposures',
                        'percentage':True,
                        'groupby':[
                                    'rcp',
                                    'rp'
                                ],
                        'climate_scenarios':['SSP1 RCP2.6','SSP5 RCP8.5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_rcp_rp.xlsx',
                        'plot_name':'percentage_by_rcp_rp.png'
                        },
                        {
                        'type':'company_type',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'rcp',
                                    'rp'
                                ],
                        'climate_scenarios':['SSP1 RCP2.6','SSP5 RCP8.5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_site_type_rcp_rp.xlsx',
                        'plot_name':'numbers_by_rcp_rp.png'
                        },
                        {
                        'type':'company_type',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'rcp',
                                    'rp'
                                ],
                        'climate_scenarios':['SSP1 RCP2.6','SSP5 RCP8.5'], 
                        'scenario_color':['#d94801','#7f2704'], 
                        'scenario_marker':['s-','^-'],     
                        'file_name':'exposures_numbers_by_site_type_rcp_rp.xlsx',
                        'plot_name':'percentage_by_rcp_rp.png'
                        },
                        {
                        'type':'company_country',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'country',
                                    'rcp',
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
                        'file_name':'exposures_numbers_by_company_country_rcp_rp.xlsx',
                        'plot_name':'max_numbers_across_hazards.png'
                        },
                        {
                        'type':'company_country',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'country',
                                    'rcp',
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
                        'file_name':'exposures_numbers_by_company_country_rcp_rp.xlsx',
                        'plot_name':'max_percentages_across_hazards.png'
                        },
                        {
                        'type':'company_name',
                        'percentage':False,
                        'groupby':[
                                    'site_type',
                                    'company_name_eng_merged',
                                    'rcp',
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
                        'file_name':'exposures_numbers_by_company_rcp_rp.xlsx',
                        'plot_name':'company_name_max_numbers_across_hazards.png'
                        },
                        {
                        'type':'company_name',
                        'percentage':True,
                        'groupby':[
                                    'site_type',
                                    'company_name_eng_merged',
                                    'rcp',
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
                        'file_name':'exposures_numbers_by_company_rcp_rp.xlsx',
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
    hazard_data_details = pd.read_csv(os.path.join(data_path,
                                    "flood_river_giri.csv"),encoding="latin1")
    hazard_keys = hazard_data_details["key"].values.tolist()
    """Get all the exposure plots
    """
    for st in range(len(plot_types)):
        if plot_types[st]["type"] in ("all_exposures","company_type"):
            exposures_df = pd.read_excel(os.path.join(exposure_results_path,
                                            f"{plot_types[st]['file_name']}"),
                                sheet_name="exposures")[
                                plot_types[st]['groupby'] + ["flood_count"]]
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

                fig, ax = plt.subplots(1,1,
                        figsize=(10,10),
                        dpi=500)
                y_min = factor*exposures["flood_count"].min()
                y_max = factor*exposures["flood_count"].max()
                rps = list(set(exposures['rp'].values)) 
                ax.plot(exposures[exposures['rcp'] == 'baseline']['rp'],
                        factor*exposures[exposures['rcp'] == 'baseline']["flood_count"],
                        'o-',color='#fd8d3c',markersize=10,linewidth=2.0,
                        label='Baseline')
                for c in range(len(plot_types[st]['climate_scenarios'])):
                    sc = plot_types[st]['climate_scenarios'][c]
                    cl = plot_types[st]['scenario_color'][c]
                    m = plot_types[st]['scenario_marker'][c]
                    exp = exposures[exposures['rcp'] == sc]
                    ax.plot(exp['rp'],
                            factor*exp["flood_count"],
                            m,color=cl,markersize=10,linewidth=2.0,
                            label=f"RCP {sc.upper()}")
                
                ax.set_xlabel('Return period (years)',fontsize=16,fontweight='bold')
                ax.set_ylabel(ylabel,fontsize=16,fontweight='bold')
                ax.set_xscale('log')
                ax.set_ylim([y_min,1.2*y_max])
                ax.tick_params(axis='both', labelsize=14)
                ax.set_xticks([t for t in rps])
                ax.set_xticklabels([str(t) for t in rps])
                ax.grid(True)         

                ax.legend(
                            loc='lower right', 
                            prop={'size':18,'weight':'bold'})
                plt.tight_layout()
                save_fig(os.path.join(figures,f"{i}_{plot_types[st]['plot_name']}"))
                plt.close()
        
        elif plot_types[st]['type'] in ("company_country","company_name"):
            exposures = pd.read_excel(os.path.join(exposure_results_path,
                                    f"{plot_types[st]['file_name']}"))
            exposures = exposures.groupby(plot_types[st]['groupby'])["flood_count"].sum().reset_index()
            site_numbers = sites.groupby(plot_types[st]['site_groupby']).size().reset_index(name='counts')
            exposures = pd.merge(exposures,site_numbers,how="left",on=plot_types[st]['site_groupby'])
            exposures["percentage_flooded"] = 100.0*exposures["flood_count"]/exposures["counts"]
            if "country" in exposures.columns.values.tolist():
                exposures["country"] = exposures["country"].str.replace("of America","")
            if plot_types[st]["percentage"] is False:
                ylabel = "Flooded number of sites"
                exposure_column = "flood_count"
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
                all_countries = [
                                    a.replace(" INDUSTRY CO.,LTD. (Tyeeli)",""
                                        ).replace(" CO., LTD",""
                                        ).replace(" Co., Ltd","") for a in all_countries
                                ]
                ac = [str(ac[:3]).upper() for ac in all_countries]
                dcsummary = [f"{a} - {b}" for (a,b) in zip(ac,all_countries)]
                dcsummary = np.reshape(dcsummary, (8, 3))
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


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
