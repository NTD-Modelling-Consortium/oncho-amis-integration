Oncho AMIS Integration
=========================

Provides the relevant wrapper code to be able to call the Oncho
simulation from the AMIS loop to compute fitted parameters.

# Usage

## Prerequistes

 * [Docker](https://docs.docker.com/desktop/)

## Setup

### Clone this repository

```shell
git clone git@github.com:NTD-Modelling-Consortium/oncho-amis-integration.git
```

### Enable SSH for Docker (skip if on macOS)
1. Add `$USER` to `docker` group

```shell
sudo usermod -a -G docker $USER

# Check that this was successful => Should print something like docker:x:111:[user]
grep docker /etc/group
```

2. Log out and log back in to change user's group ID to docker
3. Add SSH keys to agent
```shell
# This command will export important environment variables and start the SSH agent process
eval ssh-agent $SHELL

# Verify that keys have been added to the agent
ssh-add -l

# If the above command says that there are no identities available then manually add your existing private key or generate new ones and add
ssh-add path/to/private_key

ssh-add -l # This should now list the key just added
```
4. Build the Docker image
```shell
DOCKER_BUILDKIT=1 docker build --ssh default=$SSH_AUTH_SOCK . -t oncho-amis
```

## Run the Docker Container
```shell
docker run -it oncho-amis:latest <task-id>
```