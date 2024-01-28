"""China supply chain sites flood exposure stats
"""
import os
import sys
from collections import OrderedDict
import pandas as pd
import geopandas as gpd
import zipfile
import numpy as np
import ast
from summary_utils import *
from tqdm import tqdm
tqdm.pandas()

def main(config):
    processed_data_path = config['paths']['data']
    output_path = config['paths']['output']

    climate_idx_cols = ['rcp','rp'] # The climate columns

    exposure_files = ["supply_chain_sites_all"]
    asset_files = ["supply_chain_sites_all_qingnan_edited"]
    company_ids = ["company_name_eng_merged"]
    asset_ids = ["site_id"]
    indexes_columns = [["site_id","company_name_eng_merged","site_type","iso_code","country"]]
    for idx,(file,asset_file,
            company_id,asset_id,index_columns
            ) in enumerate(zip(exposure_files,asset_files,company_ids,asset_ids,indexes_columns)):
        stats_path = os.path.join(output_path,"global_exposure_stats_giri")
        if os.path.exists(stats_path) == False:
            os.mkdir(stats_path)

        hazard_data_details = pd.read_csv(os.path.join(processed_data_path,
                                        "flood_river_giri.csv"),encoding="latin1")
        hazard_keys = hazard_data_details["key"].values.tolist()

        # assets = gpd.read_file(os.path.join(processed_data_path,
        #                                 "supply_chain_study",
        #                                 f"{asset_file}.gpkg"),
        #                                 layer="nodes")
        assets = pd.read_csv(os.path.join(processed_data_path,
                                        "supply_chain_study",
                                        f"{asset_file}.csv")
                                        )
        # exposures = pd.read_parquet(os.path.join(output_path,
        #                         "sites_flood_intersections",
        #                         f"{file}_splits__aqueduct_river__nodes.geoparquet"))
        exposures = pd.read_parquet(os.path.join(output_path,
                                "sites_flood_intersections",
                                f"{file}_splits__flood_river_giri__nodes.geoparquet")) 
        exposures = exposures[[asset_id]+hazard_keys]
        exposures = pd.merge(exposures,assets[index_columns],how="left",on=[asset_id])
        exposure_stats = add_rows_and_transpose(exposures,hazard_data_details,
                                    index_columns,
                                    hazard_keys,
                                    climate_idx_cols,value_name="flood_depth_m")
        exposure_stats["flood_depth_m"] = 0.01*exposure_stats["flood_depth_m"]
        exposure_stats = exposure_stats[exposure_stats["flood_depth_m"] > 0]
        exposure_stats["flood_count"] = 1
        exposure_stats.to_csv(os.path.join(stats_path,
                            f"{file}_climate_rcp_rp_flood_depth.csv"),
                            index=False)
        exposures[hazard_keys] = np.where(exposures[hazard_keys]>0,1,0)
        exposures.to_csv(os.path.join(stats_path,
                            f"{file}_splits__flood_river_giri__nodes_flood_binary.csv"),
                            index=False)

        stats_combinations = [
                                {
                                'type':'exposure_depths',
                                'addby':[],
                                'groupby':index_columns + climate_idx_cols,
                                'file_name':'exposures_flood_depths_by_rcp_rp.xlsx',
                                'generate_quantiles':False 
                                },
                                {
                                'type':'exposure_expected_depths',
                                'addby':[],
                                'groupby':index_columns + [
                                            'rcp',
                                        ],
                                'file_name':'exposures_expected_flood_depths_by_rcp.xlsx',
                                'generate_quantiles':False 
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            'rcp',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_rcp.xlsx',
                                'generate_quantiles':False
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            'country',
                                            'rcp',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_country_rcp_rp.xlsx',
                                'generate_quantiles':False 
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            f'{company_id}',
                                            'site_type',
                                            'rcp',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_company_rcp_rp.xlsx',
                                'generate_quantiles':False 
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            'site_type',
                                            'rcp',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_site_type_rcp_rp.xlsx',
                                'generate_quantiles':False 
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            f'{company_id}',
                                            'site_type',
                                            'country',
                                            'iso_code',
                                            'rcp',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_company_country_rcp_rp.xlsx',
                                'generate_quantiles':False 
                                },
                            ]

        quantile_combinations = [
                                {
                                'type':'exposure_depths',
                                'addby':[],
                                'groupby':index_columns +[
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_flood_depths_by_rcp_rp_epoch.xlsx',
                                'generate_quantiles':True 
                                },
                                {
                                'type':'exposure_expected_depths',
                                'addby':[],
                                'groupby':index_columns +[
                                            'rcp',
                                            'epoch'
                                        ],
                                'file_name':'exposures_expected_flood_depths_by_rcp_epoch.xlsx',
                                'generate_quantiles':True 
                                },
                                {
                                'type':'exposures',
                                'groupby':[
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_rcp_rp_epoch.xlsx',
                                },
                                {
                                'type':'exposures',
                                'groupby':[
                                            'country',                                            
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_country_rcp_rp_epoch.xlsx',
                                },
                                {
                                'type':'exposures',
                                'groupby':[
                                            f'{company_id}',
                                            'site_type',
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_company_rcp_rp_epoch.xlsx',
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            'site_type',
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_site_type_rcp_rp_epoch.xlsx',
                                'generate_quantiles':True 
                                },
                                {
                                'type':'exposures',
                                'addby':[],
                                'groupby':[
                                            f'{company_id}',
                                            'site_type',
                                            'country',
                                            'iso_code',
                                            'rcp',
                                            'epoch',
                                            'rp'
                                        ],
                                'file_name':'exposures_numbers_by_company_country_rcp_rp_epoch.xlsx',
                                'generate_quantiles':True 
                                },
                            ]
        for st in range(len(stats_combinations)):
            write_values = True    

            if stats_combinations[st]['type'] != 'benefit_costs':
                output_excel = os.path.join(stats_path,
                                    f"{stats_combinations[st]['file_name']}")
                stats_wrtr = pd.ExcelWriter(output_excel)

            if stats_combinations[st]['generate_quantiles'] is True:
                output_excel = os.path.join(stats_path,
                                f"{quantile_combinations[st]['file_name']}")
                quantile_wrtr = pd.ExcelWriter(output_excel)

            if stats_combinations[st]['type'] in ('exposures','exposure_depths'):
                values_df = exposure_stats.copy()
                if stats_combinations[st]['type'] == "exposure_depths":
                    stats_by = "flood_depth_m"
                else:
                    stats_by = "flood_count"
            elif stats_combinations[st]['type'] == 'exposure_expected_depths':
                df = pd.read_csv(os.path.join(output_path,
                                "sites_flood_intersections",
                                f"{file}_expected_flood_depths.csv"))
                values_df = pd.merge(df,assets[index_columns],how="left",on=[asset_id])
                stats_by = "expected_flood_depth_undefended"
            group_and_write_results(values_df,
                stats_combinations[st]['groupby'],stats_wrtr,
                quantile_combinations[st]['groupby'],
                [stats_by],"exposures",
                write_values=write_values,
                generate_groups = True,
                generate_quantiles=stats_combinations[st]['generate_quantiles'])
            
            if stats_combinations[st]['type'] != 'benefit_costs':
                stats_wrtr.close()
            if stats_combinations[st]['generate_quantiles'] is True:
                quantile_wrtr.close()

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
