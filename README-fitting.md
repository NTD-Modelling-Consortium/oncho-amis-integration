Scripts used for Oncho fitting and near term projections
================

- On HPC clusters, scripts that take long to run must be run through Slurm. Tables below 
show shell scripts because of this.

- Most up to date code is here to do a full rerun of all IUs: <https://github.com/NTD-Modelling-Consortium/oncho-amis-integration/tree/updateFittingCode>

- My (Raiha's) Warwick directory only has the recent runs for treatment naive + a few other IUs. Evandro's runs are saved to the cloud, but note that the batch numbers here don't necessarily align with those in my Warwick directory, and some IUs from EKs fits will have been refitted.


### Preparing histories and maps for the fitting 

| R script  (see `Maps/` directory)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| multipletimepoints_preprocess_map_and_histories.R           | run_fitting_prep.sh           |

- `multipletimepoints_preprocess_map_and_histories.R`: produces the maps and the histories from 1975-2021 and identifies batches.
- Note that there is a hard-coded variable `batch_with_no_treatment` defining the batch number which contains IUs with no MDA/VC history. 

### Setup virtual environment

- On the cluster you need to load the Python module. 

  - On the Warwick cluster:

  ```
  module load  GCCcore/12.2.0 Python/3.10.8
  ```

  - Note that EPIONCHO-IBM model is currently only available for Python 3.10 (the only option currently available on BMRC cluster is Python 3.11 hence using Warwick cluster instead)

- Clone and install EPIONCHO-IBM library <https://github.com/NTD-Modelling-Consortium/EPIONCHO-IBM/> following instructions in the Github repo.

	- Note that you probably need to generate an ssh key (<https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent>) and add it to your Github (<https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>) to be able to clone the repo
	- If not already installed, the Poetry package will also need to be installed. Follow directions in README of EPIONCHO-IBM repo
	- Code below installs master branch in the EPIONCHO-IBM repo and should create virtual environment ".venv" within the EPIONCHO-IBM directory
	
```
git clone git@github.com:NTD-Modelling-Consortium/EPIONCHO-IBM.git
cd EPIONCHO-IBM
poetry config virtualenvs.in-project true
poetry install 
```

- You will also need to install the endgame_postprocessing repository (<https://github.com/NTD-Modelling-Consortium/endgame-postprocessing>) into the EPIONCHO-IBM directory

```
cd EPIONCHO-IBM
source .venv/bin/activate
pip install git+https://github.com/NTD-Modelling-Consortium/endgame-postprocessing.git
```
	

### Running the fitting

| R script                                   | Corresponding shell script    |
|:------------------------------------------------------------|:-------------|
| oncho-endgame-multipletimepts.R            | runFitting.sh                 |
| oncho-endgame-multipletimepts-sigma0.025.R | runFitting_sigma0.025.sh      |

- Make sure `r_wrapper_endgame_fitting_multipletimepts.py` is in `EPIONCHO-IBM/`. This script is the wrapper to interact with the python model.

- After **runFitting.sh** jobs are complete, identify batches where at least 1 IU didn't reach the target ESS (using code below) and use to update batch numbers in **runFitting_sigma0.025.sh**

  ```
  grep "Some locations did not reach target ESS" outputs/log/*
  ```

- After sigma=0.025 runs complete, keep note of IUs that still had ESS < 200 so we can let Igor know not to run these.


### Preparing histories for the projections

| R script  (see `Maps/` directory)                           | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| multipletimepoints_projections_inputs.R                     | run_projections_prep.sh       |
| IUsWithInsufficientESS.R                                    | NA                            |     

- `multipletimepoints_projections_inputs.R`: produces the histories from 1975-2025, produces parameter posterior samples and assigns new batches for projections.

- In the file `multipletimepoints_projections_inputs.R`, make sure you update the object `failed_ids` to include any batches that failed altogether due to numerical issues (singular matrices).

- `IUsWithInsufficientESS.R`:  finds IUs that have ESS < 200 (after also trying higher sigma=0.025)


### Plots for the model fits (not required for the pipeline, this is just to sense check the results of the fitting)

- `plot_posteriors.R`:  saves plots to summarise the fitting results and individual trajectories for each IU. Note I transferred the fitting results (`outputs/` folder) to my local machine and ran this, not on the cluster.
- Also requires ESPEN_2021 shapefiles and output from `multipletimepoints_projections_inputs.R`  


### Running the near term projections


| Python script                                               | Corresponding shell script    |
|:------------------------------------------------------------|:------------------------------|
| wrappersimulationsmultipletimepoints.py                     | runProj.sh                    |


- Make sure `wrappersimulationsmultipletimepoints.py` is located in `EPIONCHO-IBM/` 
- Indexing for the array numbers of the projections starts from 1 and goes up to the total number of IUs (or batches if batch size is >1 IU in `multipletimepoints_projections_inputs.R`)

