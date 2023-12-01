""" Intersect assets with hazard maps
"""

import os
import json
import subprocess     

def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config

def main(config):
    # Set global paths
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    results_data_path = config['paths']['output']
    
    hazard_csv = os.path.join(processed_data_path,
                                "aqueduct_river.csv")
    features_csv = os.path.join(processed_data_path,
                                "site_layers.csv")
    run_intersections = True  # Set to True is you want to run this process
    if run_intersections is True:
        args = [
                "env", 
                "SNAIL_PROGRESS=1",
                "snail",
                "-x",
                "-vv",
                "process",
                "--features",
                f"{features_csv}",
                "--rasters",
                f"{hazard_csv}",
                "--directory",
                f"{processed_data_path}"
                ]
        print ("* Start the processing of infrastrucutre hazard raster intersections")
        print (args)
        subprocess.run(args,check=True)

        print ("* Done with the processing of infrastrucutre hazard raster intersections")



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
