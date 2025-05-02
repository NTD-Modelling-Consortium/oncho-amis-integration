# syntax=docker/dockerfile:1

# https://biohpc.cornell.edu/doc/CondaInContainer.html
FROM continuumio/miniconda3:latest

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
        build-essential cmake vim curl git unzip openssh-client

RUN conda install --yes --name base -c conda-forge \
        python=3.10 \
        r-base=4.3.2 \
        r-reticulate \
        r-truncnorm   \
        r-dplyr \
        r-magrittr \
        r-invgamma \
        r-tidyr \
        r-devtools
RUN Rscript -e "install.packages(c('AMISforInfectiousDiseases'), repos='http://cran.us.r-project.org')"
RUN conda clean -a -y

RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH=/root/.local/bin:$PATH

# Verify installations
RUN R --version && python3 --version && poetry --version

ARG ONCHO_AMIS_DIR=/ntdmc/oncho-amis-integration

# Get Oncho
# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
        ssh-keyscan github.com >> ~/.ssh/known_hosts
ADD git@github.com:NTD-Modelling-Consortium/EPIONCHO-IBM.git ${ONCHO_AMIS_DIR}/model/EPIONCHO-IBM

WORKDIR ${ONCHO_AMIS_DIR}
RUN cd model/EPIONCHO-IBM && poetry install

ARG MTP_PREPROCESS_PROJECTIONS_DIR=mtp-preprocess_projections
ARG FITTING_INPUTS_URL=https://storage.googleapis.com/ntd-data-storage/pipeline/oncho/fitting-inputs

RUN mkdir -p ${MTP_PREPROCESS_PROJECTIONS_DIR}

# Copy the scripts into the image
COPY run_fit.sh ${ONCHO_AMIS_DIR}
COPY oncho-endgame-multipletimepts.R ${ONCHO_AMIS_DIR}
COPY r_wrapper_endgame_fitting_multipletimepts.py ${ONCHO_AMIS_DIR}

ADD https://storage.googleapis.com/ntd-data-storage/pipeline/oncho/fitting-inputs/Maps.tar.gz ${ONCHO_AMIS_DIR}
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/oncho/fitting-inputs/preprocess-histories-model-output.tar.gz ${MTP_PREPROCESS_PROJECTIONS_DIR}

RUN tar --no-same-owner -xzf Maps.tar.gz -C ${ONCHO_AMIS_DIR} && \
    rm ${ONCHO_AMIS_DIR}/Maps.tar.gz && \
    tar --no-same-owner -xzf ${MTP_PREPROCESS_PROJECTIONS_DIR}/preprocess-histories-model-output.tar.gz -C ${MTP_PREPROCESS_PROJECTIONS_DIR}/ && \
    rm ${MTP_PREPROCESS_PROJECTIONS_DIR}/preprocess-histories-model-output.tar.gz

COPY mtp-preprocess_projections/multipletimepoints_preprocess_histories.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/multipletimepoints_preprocess_map.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/multipletimepoints_projections_inputs.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/run_projections_inputs.sh ${MTP_PREPROCESS_PROJECTIONS_DIR}/

ENV ONCHO_AMIS_DIR=${ONCHO_AMIS_DIR}

ENTRYPOINT ["bash", "run_fit.sh"]