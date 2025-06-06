# syntax=docker/dockerfile:1

FROM condaforge/miniforge3

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

ARG ONCHO_AMIS_DIR=/ntdmc/oncho-amis-integration
ARG ONCHO_MODEL_DIR=${ONCHO_AMIS_DIR}/model/EPIONCHO-IBM
ARG FITTING_PREP_DIR=${ONCHO_AMIS_DIR}/fitting-prep
ARG FITTING_DIR=${ONCHO_AMIS_DIR}/fitting
ARG PROJECTIONS_PREP_DIR=${ONCHO_AMIS_DIR}/projections-prep
ARG PROJECTIONS_DIR=${ONCHO_AMIS_DIR}/projections

RUN apt update && apt install -y \
	build-essential \
	cmake \
	vim \
	curl \
	git \
	unzip \
	openssh-client

RUN conda install --override-channels -c conda-forge -c r --yes --name base \
	python=3.10 \
	poetry \
	r-base \
	r-reticulate \
	r-truncnorm   \
	r-dplyr \
	r-magrittr \
	r-invgamma \
	r-tidyr \
	r-devtools \
	r-openxlsx \
	r-hmisc \
	r-mclust \
	r-mnormt \
	r-optparse \
	r-rcpp \
	r-rcpparmadillo \
	r-mice

RUN Rscript -e "install.packages('weights', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library/')"
RUN Rscript -e "install.packages('AMISforInfectiousDiseases', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library/')"
RUN conda clean -a -y

# Verify installations
RUN R --version && python3 --version && poetry --version

# Get Oncho
# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
	ssh-keyscan github.com >> ~/.ssh/known_hosts
ADD git@github.com:NTD-Modelling-Consortium/EPIONCHO-IBM.git ${ONCHO_MODEL_DIR}
ADD git@github.com:NTD-Modelling-Consortium/endgame-postprocessing.git ${ONCHO_AMIS_DIR}/endgame-postprocessing

WORKDIR ${ONCHO_AMIS_DIR}
RUN cd ${ONCHO_MODEL_DIR} && \
	poetry install && \
	poetry run pip install -e ${ONCHO_AMIS_DIR}/endgame-postprocessing

RUN mkdir -p ${FITTING_PREP_DIR}/{inputs,artefacts,scripts} && \
	mkdir -p ${FITTING_PREP_DIR}/artefacts/{Maps,model_output} && \
	mkdir -p ${FITTING_DIR}/{inputs,artefacts,scripts} && \
	mkdir -p ${PROJECTIONS_PREP_DIR}/{inputs,artefacts,scripts} && \
	mkdir -p ${PROJECTIONS_DIR}/{inputs,artefacts,scripts}

# Copy the scripts into the container
# Order of the shell scripts indicates the stages of the pipeline
# ⬇ Stage 1: Fitting-Prep
COPY mtp-preprocess_projections/multipletimepoints_preprocess_map_and_histories.R ${FITTING_PREP_DIR}/scripts/
COPY run_fitting_prep.sh ${FITTING_PREP_DIR}/scripts/
# ⬇ Stage 2: Fitting
COPY r_wrapper_endgame_fitting_multipletimepts.py ${FITTING_DIR}/scripts/
COPY oncho-endgame-multipletimepts.R ${FITTING_DIR}/scripts/
COPY run_fitting.sh ${FITTING_DIR}/scripts/
# ⬇ Stage 3: Projections-Prep
COPY mtp-preprocess_projections/multipletimepoints_projections_inputs.R ${PROJECTIONS_PREP_DIR}/scripts/
COPY run_projections_inputs.sh ${PROJECTIONS_PREP_DIR}/scripts/
COPY ius_with_insufficient_ess.r ${PROJECTIONS_PREP_DIR}/scripts/
# ⬇ Stage 4: Projections
COPY wrapper-simulations-multipletimepoints.py ${PROJECTIONS_DIR}/scripts/
COPY run_projections_to_2026.sh ${PROJECTIONS_DIR}/scripts/

COPY run_pipeline.py ${ONCHO_AMIS_DIR}/

# Download the artefacts created by the fitting prep script
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/oncho/fitting-inputs/artefacts-fitting-prep.tar.gz ${ONCHO_AMIS_DIR}
RUN tar --no-same-owner -xzf artefacts-fitting-prep.tar.gz -C ${ONCHO_AMIS_DIR} && \
	mv ${ONCHO_AMIS_DIR}/artefacts-fitting-prep/Maps/{maps_joint,ALL_prevalence_map_multipletimespoints.rds,iu_task_lookup.rds} ${FITTING_PREP_DIR}/artefacts/Maps && \
	mv ${ONCHO_AMIS_DIR}/artefacts-fitting-prep/model_output/fitting/* ${FITTING_PREP_DIR}/artefacts/model_output && \
	mv ${ONCHO_AMIS_DIR}/artefacts-fitting-prep/Maps/{ALL_prevalence_map.rds,Full_histories_df_popinfo_ALL_minimal_070425_listlabels.rds} ${FITTING_PREP_DIR}/inputs && \
	rm ${ONCHO_AMIS_DIR}/artefacts-fitting-prep.tar.gz && \
	rm -rf ${ONCHO_AMIS_DIR}/artefacts-fitting-prep

ENV ONCHO_AMIS_DIR=${ONCHO_AMIS_DIR}
ENV ONCHO_MODEL_DIR=${ONCHO_MODEL_DIR}
ENV PATH_TO_FITTING_PREP_ARTEFACTS="$ONCHO_AMIS_DIR/fitting-prep/artefacts"
ENV PATH_TO_FITTING_PREP_INPUTS="$ONCHO_AMIS_DIR/fitting-prep/inputs"
ENV PATH_TO_FITTING_PREP_SCRIPTS="$ONCHO_AMIS_DIR/fitting-prep/scripts"

ENV PATH_TO_FITTING_ARTEFACTS="$ONCHO_AMIS_DIR/fitting/artefacts"
ENV PATH_TO_FITTING_INPUTS="$ONCHO_AMIS_DIR/fitting/inputs"
ENV PATH_TO_FITTING_SCRIPTS="$ONCHO_AMIS_DIR/fitting/scripts"

ENV PATH_TO_PROJECTIONS_PREP_ARTEFACTS="$ONCHO_AMIS_DIR/projections-prep/artefacts"
ENV PATH_TO_PROJECTIONS_PREP_INPUTS="$ONCHO_AMIS_DIR/projections-prep/inputs"
ENV PATH_TO_PROJECTIONS_PREP_SCRIPTS="$ONCHO_AMIS_DIR/projections-prep/scripts"

ENV PATH_TO_PROJECTIONS_ARTEFACTS="$ONCHO_AMIS_DIR/projections/artefacts"
ENV PATH_TO_PROJECTIONS_INPUTS="$ONCHO_AMIS_DIR/projections/inputs"
ENV PATH_TO_PROJECTIONS_SCRIPTS="$ONCHO_AMIS_DIR/projections/scripts"

ENV PATH_TO_MAPS="$PATH_TO_FITTING_PREP_ARTEFACTS/Maps"
ENV PATH_TO_FITTING_PREP_MODEL_OUTPUT="$PATH_TO_FITTING_PREP_ARTEFACTS/model_output"
ENV PATH_TO_PROJECTIONS_PREP_MODEL_OUTPUT="$PATH_TO_PROJECTIONS_PREP_ARTEFACTS/model_output"
ENV PATH_TO_OUTPUTS="$PATH_TO_FITTING_ARTEFACTS/outputs"
ENV PATH_TO_MTP_PREPROCESS_PROJ_SCRIPTS="$ONCHO_AMIS_DIR/mtp-preprocess_projections"
ENV PATH_TO_MODEL="$ONCHO_MODEL_DIR"

ENTRYPOINT ["python", "run_pipeline.py"]