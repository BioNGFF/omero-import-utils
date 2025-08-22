#!/bin/bash

# Usage: metadata/test_register.sh <project_id>

# This script creates a new Dataset and imports various Zarr images into it.

# NB: There are also 2 plate examples at the end of this script, commented out for now.

if [ $# -ne 1 ]; then
  echo "Usage: $0 <project_id>"
  exit 1
fi

P="$1"

# create a Dataset in OMERO. Named by Timestamp, linked to Project
name=$(date "+%Y:%m:%d_%H:%M:%S")
dataset=$(omero obj new Dataset name="$name")
omero obj new ProjectDatasetLink parent=Project:$P child=$dataset

echo $dataset

# Classic Sample - duplicated below
python register.py --target "$dataset" https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0062A/6001240.zarr

# Embassy image - duplicate:
python register.py --target "$dataset" https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0101A/13457227.zarr
python register.py --target "$dataset" --endpoint https://uk1s3.embassy.ebi.ac.uk/ --nosignrequest s3://idr/zarr/v0.4/idr0101A/13457227.zarr

# AWS images - first 2 are duplicates:
python register.py --target "$dataset" https://s3.us-east-1.amazonaws.com/gs-public-zarr-archive/CMU-1.ome.zarr/0
python register.py --target "$dataset" --endpoint https://s3.us-east-1.amazonaws.com/ s3://gs-public-zarr-archive/CMU-1.ome.zarr/0 --nosignrequest
# data.source - flaky!
python register.py --target "$dataset" https://data.source.coop/joshmoore/idr-ome-ngff-samples/v0.4/idr0062A/6001240.zarr

# Google cloud - bioformats2raw.layout - 3 images
python register.py --target "$dataset" https://storage.googleapis.com/jax-public-ngff/example_v2/LacZ_ctrl.zarr

# bioformats2raw.layout (single image)
python register.py --target "$dataset" https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0048A/9846151.zarr/

# Test rendering settings (image also has labels)
python register.py --target "$dataset" https://uk1s3.embassy.ebi.ac.uk/bia-idr-integration/S-BIAD1961/experimentA/RP-TMA/RP-TMA-3_ROI01_multichannel.zarr/experimentA_backup/RP-TMA/RP-TMA-3_ROI01_multichannel.zarr/

# plate from idr0011 (48 Wells, 1 field):
# python register.py https://uk1s3.embassy.ebi.ac.uk/bia-integrator-data/S-BIAD866/7f95aba3-cfbf-4ae8-a106-edaa36f5b07f/7f95aba3-cfbf-4ae8-a106-edaa36f5b07f.zarr/

# plate from idr0001 (96 Wells, 6 fields - 6 acquisitions):
# python register.py https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0001A/2551.zarr
