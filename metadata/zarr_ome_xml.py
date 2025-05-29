# requires Zarr v2 and ome-types
# pip install ome-types zarr==2.18.7

from zarr.storage import FSStore
import json
import dask.array as da
import sys
from ome_types import to_xml
from ome_types.model import OME, Image, Pixels, Plate, Well, WellSample, ImageRef
from ome_types.model.simple_types import ImageID, PixelsID, PixelType, PlateID, WellID, WellSampleID


pixels_id_counter = 1
image_id_counter = 1
well_id_counter = 1
well_sample_id_counter = 1


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
            array_path =  f"{well_sample_path}/{image_json["multiscales"][0]["datasets"][0]["path"]}"
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


def handle_image(store, zarr_uri):
    zattrs = json.loads(store.get(".zattrs"))
    array_path = zattrs["multiscales"][0]["datasets"][0]["path"]
    array_data = da.from_zarr(store, array_path)

    sizes = {}
    shape = array_data.shape
    axes = zattrs["multiscales"][0]["axes"]
    for dim, size in zip(axes, shape):
        sizes[dim["name"]] = size
    pixels_type = array_data.dtype.name

    object_name = zarr_uri.rstrip("/").split("/")[-1].split(".")[0]

    img = create_image(object_name, pixels_type, sizes)
    ome = OME()
    ome.images.append(img)
    write_xml(ome, object_name)    


if len(sys.argv) < 2:
    print("Error: Please provide the Zarr URI as a command line argument")
    print("Usage: python zarr_ome_xml.py <zarr_uri>")
    sys.exit(1)

zarr_uri = sys.argv[1]
zattrs = json.load(open(f"{zarr_uri}/.zattrs"))

if "plate" in zattrs:
    store = FSStore(zarr_uri)
    handle_plate(store, zarr_uri)
else:
    if "bioformats2raw.layout" in zattrs and zattrs["bioformats2raw.layout"] == 3:
        store = FSStore(f"{zarr_uri}/0")
    else:
        store = FSStore(zarr_uri)
    handle_image(store, zarr_uri)
