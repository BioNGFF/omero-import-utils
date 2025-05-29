#!/bin/bash

# Can be called like "omero import ... SOME_IMAGE", ie. additional arguments will be
# passed through.
# E.g. "./import.sh -T Dataset:18651 /data/testing/cat.ome.zarr/OME/METADATA.ome.xml"

# Assumes 'omero' micromamba environment (with omero-cli-render and omero-cli-zarr)
# is installed.

# Activate OMERO environment
eval "$(micromamba shell hook -s bash)"
micromamba activate omero

# Check if reader.txt exists and create it if it doesn't
if [ ! -f "/tmp/reader.txt" ]; then
    echo "loci.formats.in.OMEXMLReader" > /tmp/reader.txt
fi

# Extract image name from arguments
for arg in "$@"; do
    if [[ "$arg" == *".ome.zarr"* ]]; then
        IMAGE_NAME="$(basename "${arg%/OME*}")"
        break
    fi
done

# Run import and process each imported image
omero import -l /tmp/reader.txt -n "${IMAGE_NAME}" $@ | grep -e "Image:[0-9]\+" -e "Plate:[0-9]\+" | while read -r image; do
    omero zarr extinfo --set "$image"
    omero render test --thumb "$image"
done
