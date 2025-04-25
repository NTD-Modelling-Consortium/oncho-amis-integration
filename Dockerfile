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
RUN Rscript -e 'devtools::install_github("drsimonspencer/AMISforInfectiousDiseases")'
RUN conda clean -a -y

RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH=/root/.local/bin:$PATH

# Verify installations
RUN R --version
RUN python3 --version
RUN poetry --version

# Get Oncho
RUN mkdir -p -m 0600 ~/.ssh && \
        ssh-keyscan github.com >> ~/.ssh/known_hosts
ADD git@github.com:NTD-Modelling-Consortium/EPIONCHO-IBM.git /ntdmc/EPIONCHO-IBM

WORKDIR /ntdmc
RUN cd EPIONCHO-IBM && poetry install

ARG ONCHO_AMIS_DIR=/ntdmc/oncho-amis-integration
ARG MTP_PREPROCESS_PROJECTIONS_DIR=${ONCHO_AMIS_DIR}/mtp-preprocess_projections

RUN mkdir -p ${ONCHO_AMIS_DIR} ${MTP_PREPROCESS_PROJECTIONS_DIR}

# Copy the scripts into the image
COPY run_fit.sh ${ONCHO_AMIS_DIR}
COPY oncho-endgame-multipletimepts.R ${ONCHO_AMIS_DIR}
COPY r_wrapper_endgame_fitting_multipletimepts.py ${ONCHO_AMIS_DIR}
COPY Maps ${ONCHO_AMIS_DIR}/Maps
COPY mtp-preprocess_projections/multipletimepoints_preprocess_histories.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/multipletimepoints_preprocess_map.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/multipletimepoints_projections_inputs.R ${MTP_PREPROCESS_PROJECTIONS_DIR}/
COPY mtp-preprocess_projections/run_projections_inputs.sh ${MTP_PREPROCESS_PROJECTIONS_DIR}/
ADD preprocess-histories-model-output.tar.gz ${MTP_PREPROCESS_PROJECTIONS_DIR}/

# ARG ONCHO_PYVENV=/root/.cache/pypoetry/virtualenvs/epioncho-ibm-U5jSXRn5-py3.10
# RUN --mount=type=cache,target=/root/.cache/R source ${ONCHO_PYVENV}/bin/activate
WORKDIR ${ONCHO_AMIS_DIR}


ENTRYPOINT ["bash"]