#!/bin/bash

# Note that due to difficulties with the private GitLab Container Registry,
# we've decided to use Docker Hub, so this script should no longer be used.
# Instead, refer to https://docs.docker.com/docker-hub/, and make sure to tag
# the image using a tag of the form vX. See the Docker Hub page at
# https://hub.docker.com/r/jcking1034/sp-server-cicd

# A small script to build the image described by 'gitlab-ci.Dockerfile', which
# can be used in '.gitlab-ci.yml'

FILE_NAME="gitlab-ci.Dockerfile"

REGISTRY_LOCATION="git.elfin.ucla:5050"
IMAGE_PATH="${REGISTRY_LOCATION}/science-processing/sp-server"

DATE_STR=$(date -u +"%y%m%d_%H%M%S")
IMAGE_NAME="gitlab-image"

FULL_DATE_STR_NAME="${IMAGE_PATH}/${IMAGE_NAME}:${DATE_STR}"
FULL_LATEST_NAME="${IMAGE_PATH}/${IMAGE_NAME}:latest"

# https://stackoverflow.com/a/36398336
docker build \
    --tag "$FULL_DATE_STR_NAME" \
    --tag "$FULL_LATEST_NAME" \
    --file "${FILE_NAME}" \
    .

# https://docs.docker.com/registry/insecure/
# https://stackoverflow.com/q/42211380
echo -e "\n\n------------------------------------------"
echo -e "To push this image to the GitLab registry:"
echo -e "- Ensure your docker daemon is updated"
echo -e "- Ensure your docker daemon has \"${REGISTRY_LOCATION}\" as an insecure registry"
echo -e "- Log in to the Gitlab registry:\n\t\$ docker login ${REGISTRY_LOCATION}"
echo -e "- Push the image tagged with the date string:\n\t\$ docker push ${FULL_DATE_STR_NAME}"
echo -e "- Push the latest image:\n\t\$ docker push ${FULL_LATEST_NAME}"
echo -e "\nIf using for GitLab CI, make sure to update 'image' in .gitlab-ci.yml!"
