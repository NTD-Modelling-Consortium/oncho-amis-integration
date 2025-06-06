# Oncho Fitting and Near-Term Projections (Historical Simulations)
## Instructions for running in the Cloud using Docker

#### Things to note:
- Raiha's Warwick directory only has the recent runs for treatment naive + a few other IUs. Evandro's runs are saved to the cloud, but note that the batch numbers here don't necessarily align with those in my Warwick directory, and some IUs from EKs fits will have been refitted.

### Setup
- Clone this repo - [oncho-amis-integration](https://github.com/NTD-Modelling-Consortium/oncho-amis-integration.git).
- At the root of the repo, build the Docker image:
    ```shell
    DOCKER_BUILDKIT=1 docker build --ssh default=$SSH_AUTH_SOCK . -t oncho-amis-pipeline
    ```
    This assumes that you have added the SSH keys on your system to your Github account.

### Pipeline
The pipeline has four stages, each corresponding to some scripts - 
1. **Fitting Preparation**
    - `multipletimepoints_preprocess_map_and_histories.R` - produces the maps and histories from 1975-2021 and identifies batches.
    - `run_fitting_prep.sh`
    
    **NOTE**: There is a hard-coded variable `batch_with_no_treatment` defining the batch number which contains IUs with no MDA/VC history.
2. **Fitting**
    - `oncho-endgame-multipletimepts.R`
    - `r_wrapper_endgame_fitting_multipletimepts.py`
    - `run_fitting.sh`

    **NOTE**:
    - At the end of the `fitting` stage, identify batches where at least one IU didn't reach the target ESS (see `fitting/artefacts/summary.csv`) and use them to run the fitting stage again but with `amis-sigma=0.025`.
    - If there are batches with IUs where `ESS < 200`, even after using `amis-sigma=0.025` then note them down so that Igor doesn't run them.
3. **Projections Preparation**
    - `multipletimepoints_projections_inputs.R` - produces the histories from 1975-2025, produces parameter samples and assigns new batches for projections.
    - `run_projections_inputs.sh`
    - `ius_with_insufficient_ess.r` - finds IUs that have `ESS < 200` (after also trying a higher `amis-sigma` value)

    **NOTE**:
    - At the end of the stage, identify batches where at least one IU didn't reach the target ESS and use them to rerun, but with `amis-sigma=0.025`.
    - If there are batches with IUs where `ESS < 200`, even after using `amis-sigma=0.025` then note them down so that Igor doesn't run them.
4. **Projections**
    - `wrapper-simulations-multipletimepoints.py`
    - `run_projections_to_2026.sh`

#### Usage
**NOTE**:
- It is expected that the shell scripts in each stage will be run, since they setup the Python virtual environment correctly.

The shell script `run_container.sh` is provided, as a wrapper around `run_pipeline.py`, for convenience to run the pipeline end-to-end, and through all the stages. In normal usage, this is the only script that is needed.

It's usage is as follows - 
```shell
usage: run_pipeline.py [-h] -i ID [--failed-ids FAILED_IDS]
                       [--stage {fitting-prep,fitting,projections-prep,nearterm-projections,all,skip-fitting-prep}]
                       [--amis-sigma AMIS_SIGMA]
                       [--amis-target-ess AMIS_TARGET_ESS]
                       [--amis-n-samples AMIS_N_SAMPLES]
                       [--amis-max-iters AMIS_MAX_ITERS]
                       [--ess-threshold ESS_THRESHOLD]

Run the oncho AMIS pipeline end-to-end

options:
  -h, --help            show this help message and exit
  -i ID, --id ID        Batch/task ID to process
  --failed-ids          FAILED_IDS
                        Comma-separated list ('id1,id2,id3...') of failed
                        batch/task IDs to skip. Only used when --id is not
                        specified and in the Projections-Prep stage.
  --stage               {fitting-prep,
                         fitting,
                         projections-prep,
                         nearterm-projections,
                         all,skip-fitting-prep}
                        Stage of the pipeline to run.
                        Options: 'fitting-prep',
                                 'fitting',
                                 'projections-prep',
                                 'nearterm-projections',
                                 'all',
                                 'skip-fitting-prep'.
                                Default is 'skip-fitting-prep'.
                        If 'all' is specified, it runs all stages in order.
                        If 'skip-fitting-prep' is specified, it skips
                        the fitting-prep stage.
  --amis-sigma          AMIS_SIGMA
                        AMIS 'sigma' parameter (default: 0.0025)
  --amis-target-ess     AMIS_TARGET_ESS
                        Target ESS parameter for AMIS (default: 500)
  --amis-n-samples      AMIS_N_SAMPLES
                        Number of AMIS samples (default: 500)
  --amis-max-iters      AMIS_MAX_ITERS
                        Maximum number of AMIS iterations (default: 50)
  --ess-threshold       ESS_THRESHOLD
                        ESS threshold parameter (default: 200)
```

Example (Default Parameters),
```shell
docker run oncho-amis-pipeline:latest --stage=skip-fitting-prep --id=11
```

This will produce the artefacts from each stage and copy them over to the host in the `artefacts` directory.
