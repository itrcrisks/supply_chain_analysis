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

    climate_idx_cols = ['model','rcp', 'epoch','rp'] # The climate columns

    exposure_files = ["polygon_sites_all"]
    asset_files = ["polygon_sites"]
    company_ids = ["name_eng"]
    asset_ids = ["site_id"]
    indexes_columns = [["site_id","name_eng","flooded_area_sqm"]]
    for idx,(file,asset_file,
            company_id,asset_id,index_columns
            ) in enumerate(zip(exposure_files,asset_files,company_ids,asset_ids,indexes_columns)):
        stats_path = os.path.join(output_path,"global_exposure_stats")
        if os.path.exists(stats_path) == False:
            os.mkdir(stats_path)

        hazard_data_details = pd.read_csv(os.path.join(processed_data_path,
                                        "aqueduct_river.csv"),encoding="latin1")
        hazard_keys = hazard_data_details["key"].values.tolist()

        assets = gpd.read_file(os.path.join(processed_data_path,
                                        "supply_chain_study",
                                        f"{asset_file}.gpkg"),
                                        layer="areas")
        assets = assets.to_crs(epsg=3995)
        assets["total_area_sqm"] = assets.geometry.area
        exposures = gpd.read_parquet(os.path.join(output_path,
                                "sites_flood_intersections",
                                f"{file}_splits__aqueduct_river__areas.geoparquet"))
        exposures = exposures.set_crs(epsg=4326)
        exposures = exposures.to_crs(epsg=3995)
        exposures["flooded_area_sqm"] = exposures.geometry.area
        exposures.drop("geometry",axis=1,inplace=True)
        exposures.to_csv(os.path.join(stats_path,
                            f"{file}_splits__aqueduct_river__areas_flood_depths.csv"),index=False) 

        exposures = exposures[index_columns+hazard_keys]
        # exposures = pd.merge(exposures,assets[index_columns],how="left",on=[asset_id])
        exposure_stats = add_rows_and_transpose(exposures,hazard_data_details,
                                    index_columns,
                                    hazard_keys,
                                    climate_idx_cols,value_name="flood_depth_m")
        exposure_stats = exposure_stats[exposure_stats["flood_depth_m"] > 0]
        exposure_stats = exposure_stats.groupby(
                            ["site_id","name_eng"] + climate_idx_cols + ["flood_depth_m"]
                            )["flooded_area_sqm"].sum().reset_index()
        exposure_stats = pd.merge(exposure_stats,
                        assets[[asset_id,"total_area_sqm"]],
                        how="left",on=[asset_id])
        exposure_stats["percentage_area_flooded"] = 100.0*exposure_stats["flooded_area_sqm"]/exposure_stats["total_area_sqm"]
        exposure_stats["flood_count"] = 1
        exposure_stats.to_csv(os.path.join(stats_path,
                            f"{file}_climate_model_rcp_rp_epoch_flood_depth.csv"),
                            index=False)
        exposures[hazard_keys] = np.where(exposures[hazard_keys]>0,1,0)
        exposures.to_csv(os.path.join(stats_path,
                            f"{file}_splits__aqueduct_river__nodes_flood_binary.csv"),
                            index=False)
        exposure_areas = exposure_stats.groupby(
                            ["site_id","name_eng"] + climate_idx_cols
                            )[["flooded_area_sqm","percentage_area_flooded"]].sum().reset_index()
        exposure_areas.to_csv(os.path.join(stats_path,
                            f"{file}_climate_model_rcp_rp_epoch.csv"),
                            index=False)
        exposure_areas = quantiles(exposure_areas,
                                ["site_id","name_eng","rcp","epoch","rp"],
                                ["flooded_area_sqm","percentage_area_flooded"])
        exposure_areas.to_csv(os.path.join(stats_path,
                            f"{file}_climate_rcp_rp_epoch.csv"),
                            index=False)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
