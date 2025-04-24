# syntax=docker/dockerfile:1

# https://biohpc.cornell.edu/doc/CondaInContainer.html
FROM continuumio/miniconda3:latest

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
        build-essential cmake vim curl git unzip openssh-client

RUN conda install --yes --name base -c conda-forge \
        python=3.10 \
        r-base=4.3.2
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
RUN mkdir -p /ntdmc/
RUN --mount=type=ssh git clone git@github.com:NTD-Modelling-Consortium/EPIONCHO-IBM.git /ntdmc/EPIONCHO-IBM

WORKDIR /ntdmc

ENTRYPOINT [ "bash" ]