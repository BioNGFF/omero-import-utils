# requires Zarr v2 and ome-types
# pip install ome-types zarr==2.18.7

import zarr
from zarr.storage import FSStore
import json
import dask.array as da
import os
import sys
from ome_types import to_xml
from ome_types.model import OME, Image, Pixels, Plate, Well, WellSample, ImageRef
from ome_types.model.simple_types import ImageID, PixelsID, PixelType, PlateID, WellID, WellSampleID

"""
Usage: python zarr_ome_xml.py <zarr_uri>

If the zarr_uri points to a plate, we create a Plate with Wells and WellSamples.

If the zarr_uri points to a single image (multiscales in /.zattrs) we create
an OME XML file with a single Image.
Otherwise, we look for zarr_uri/0, zarr_uri/1 etc and create an OME XML file
with all the images found, at zarr_uri/OME/METADATA.ome.xml
"""

pixels_id_counter = 0
image_id_counter = 0
well_id_counter = 0
well_sample_id_counter = 0


def write_xml(ome, name):
    xml_text = to_xml(ome)
    with open(f'{name}.ome.xml', 'w') as f:
        f.write(xml_text)


def create_image(image_name, pixels_type, sizes):
    global pixels_id_counter
    global image_id_counter
    pixels_id = pixels_id_counter
    pixels_id_counter += 1
    image_id = image_id_counter
    image_id_counter += 1

    ptype = PixelType(pixels_type)
    pixels = Pixels(
        id=PixelsID("Pixels:%s" % pixels_id),
        dimension_order="XYZCT",
        size_c=sizes.get("c", 1),
        size_t=sizes.get("t", 1),
        size_z=sizes.get("z", 1),
        size_x=sizes.get("x", 1),
        size_y=sizes.get("y", 1),
        type=ptype,
        metadata_only=True,
    )
    return Image(id=ImageID("Image:%s" % image_id), pixels=pixels, name=image_name)
    

def handle_plate(store, zarr_uri):
    global well_id_counter
    global well_sample_id_counter
    zattrs = json.loads(store.get(".zattrs"))
    plate_attrs = zattrs["plate"]
    n_cols = len(plate_attrs["columns"])
    n_rows = len(plate_attrs["rows"])
    object_name = zarr_uri.rstrip("/").split("/")[-1].split(".")[0]

    ome = OME()
    plate = Plate(id=PlateID("Plate:1"), name=object_name, columns=n_cols, rows=n_rows)
    
    for well_attrs in plate_attrs["wells"]:
        row_index = well_attrs["rowIndex"]
        col_index = well_attrs["columnIndex"]
        well = Well(id=WellID("Well:%s" % well_id_counter), row=row_index, column=col_index)
        well_id_counter += 1
        
        well_path = well_attrs['path']
        well_samples_json = json.loads(store.get(f"{well_path}/.zattrs"))
        for index, sample_attrs in enumerate(well_samples_json["well"]["images"]):
            well_sample = WellSample(id=WellSampleID("WellSample:%s" % well_sample_id_counter), index=index)
            well_sample_id_counter += 1
            well_sample_path = f"{well_path}/{sample_attrs['path']}"
            image_json =  json.loads(store.get(f"{well_sample_path}/.zattrs"))
            array_path =  f"{well_sample_path}/{image_json['multiscales'][0]['datasets'][0]['path']}"
            array_data = da.from_zarr(store, array_path)
            sizes = {}
            shape = array_data.shape
            axes = image_json["multiscales"][0]["axes"]
            for dim, size in zip(axes, shape):
                sizes[dim["name"]] = size
            pixels_type = array_data.dtype.name
            img = create_image(well_sample_path, pixels_type, sizes)
            ome.images.append(img)
            well_sample.image_ref = ImageRef(id=img.id)
            well.well_samples.append(well_sample)

        plate.wells.append(well)

    ome.plates.append(plate)
    write_xml(ome, object_name) 


def handle_images(zarr_uri, zattrs):

    ome = OME()

    if zattrs is not None and "multiscales" in zattrs:
        image_uris = [zarr_uri]
    else:
        # look for /0, /1 etc
        image_uris = []
        index = 0
        while True:
            image_uri = os.path.join(zarr_uri, str(index))
            print("Checking for", image_uri)
            if not os.path.exists(image_uri):
                break
            image_uris.append(image_uri)
            index += 1
    if len(image_uris) == 0:
        print("No images found in", zarr_uri)

    print("image_uris", image_uris)
    for image_uri in image_uris:
        print("Handling", image_uri)
        store = FSStore(image_uri)
        zattrs = json.loads(store.get(".zattrs"))
        array_path = zattrs["multiscales"][0]["datasets"][0]["path"]
        array_data = da.from_zarr(store, array_path)
        object_name = image_uri.rstrip("/").split("/")[-1].split(".")[0]
        if zattrs["multiscales"][0].get("name", "/") != "/":
            object_name = zattrs["multiscales"][0]["name"]

        sizes = {}
        shape = array_data.shape
        axes = zattrs["multiscales"][0]["axes"]
        for dim, size in zip(axes, shape):
            sizes[dim["name"]] = size
        pixels_type = array_data.dtype.name

        img = create_image(object_name, pixels_type, sizes)
        ome.images.append(img)

    out_path = os.path.join(zarr_uri, "OME", "METADATA")
    os.makedirs(os.path.join(zarr_uri, "OME"))
    write_xml(ome, out_path)

    # If the zarr_url directory isn't a zarr group, we create it
    if not os.path.exists(os.path.join(zarr_uri, ".zgroup")):
        zgroup = zarr.open_group(zarr_uri, mode='a')
        zgroup.attrs["bioformats2raw.layout"] = 3


if len(sys.argv) < 2:
    print("Error: Please provide the Zarr URI as a command line argument")
    print("Usage: python zarr_ome_xml.py <zarr_uri>")
    sys.exit(1)

print("ARGS", sys.argv)

zarr_uri = sys.argv[1]

zattrs = None
plate_uri = None
try:
    zattrs = json.load(open(f"{zarr_uri}/.zattrs"))
    if "plate" in zattrs:
        plate_uri = zarr_uri
except Exception as e:
    print("Not a plate:", e)
    pass

if plate_uri is not None:
    print("Handling plate", plate_uri)
    store = FSStore(plate_uri)
    handle_plate(store, plate_uri)
else:
    # If zattrs not found, we look in /0 etc for images
    handle_images(zarr_uri, zattrs)
