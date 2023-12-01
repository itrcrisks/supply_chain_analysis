# Supply Chain Analysis

Scripts to do intersections of point and polygon datasets of global sites with river flood maps.

The flood datasets are obtained from [WRI Aqueduct](https://www.wri.org/data/aqueduct-floods-hazard-maps)
They can be downloaded by following the instructions [here](https://github.com/nismod/aqueduct)

## Development setup

Clone this repository:

    git clone https://github.com/itrcrisks/supply_chain_analysis.git

Move into the cloned folder:

    cd supply_chain_analysis

Create a conda environment using
[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
to install packages specified in the [`environment.yaml`
file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually):

    micromamba env create -f environment.yaml

Activate it:

    micromamba activate supply_chain_analysis

Set-up the data folder paths in the `conig.json` file, dervied from the `conig.template.json` 


