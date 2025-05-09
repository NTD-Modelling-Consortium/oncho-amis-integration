#!/bin/bash

# This script runs the Docker container for the Oncho AMIS integration project.
# It mounts the necessary directories for the Maps and model output.
# Make sure to run this script from the root directory of the project.
# Usage: ./run_docker_container.sh
# Check if Docker is installed
if ! command -v docker &> /dev/null
then
    echo "Docker could not be found. Please install Docker to run this script."
    exit
fi
# Check if Docker is running
if ! docker info &> /dev/null
then
    echo "Docker is not running. Please start Docker to run this script."
    exit
fi
# Check if the Docker image is built
if ! docker image inspect oncho-amis &> /dev/null
then
    echo "Docker image 'oncho-amis' not found. Please build the Docker image before running this script."
    exit
fi

docker volume create oncho-amis-volume

# Run the Docker container
docker run --entrypoint bash \
        --rm \
        -v "$(pwd)/Maps:/tmp/Maps" \
        -v "$(pwd)/mtp-preprocess_projections:/tmp/model_output/" \
        -v "$(pwd)/outputs:/tmp/outputs" \
        -it oncho-amis:latest