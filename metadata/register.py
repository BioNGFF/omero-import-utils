#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2025 University of Dundee & Open Microscopy Environment.
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from urllib.parse import urlsplit
import numpy as np
import zarr

import argparse

from numpy import iinfo, finfo

# from getpass import getpass

from omero.cli import cli_login
from omero.gateway import BlitzGateway, ColorHolder
from omero.gateway import OMERO_NUMPY_TYPES

import omero
import omero_rois

from omero.model.enums import PixelsTypeint8, PixelsTypeuint8, PixelsTypeint16
from omero.model.enums import PixelsTypeuint16, PixelsTypeint32
from omero.model.enums import PixelsTypeuint32, PixelsTypefloat
from omero.model.enums import PixelsTypecomplex, PixelsTypedouble

from omero.model import ExternalInfoI, RoiI, MaskI
from omero.rtypes import rbool, rdouble, rint, rlong, rstring, rtime


AWS_DEFAULT_ENDPOINT = "s3.us-east-1.amazonaws.com"

PIXELS_TYPE = {'int8': PixelsTypeint8,
               'int16': PixelsTypeint16,
               'uint8': PixelsTypeuint8,
               'uint16': PixelsTypeuint16,
               'int32': PixelsTypeint32,
               'uint32': PixelsTypeuint32,
               'float_': PixelsTypefloat,
               'float8': PixelsTypefloat,
               'float16': PixelsTypefloat,
               'float32': PixelsTypefloat,
               'float64': PixelsTypedouble,
               'complex_': PixelsTypecomplex,
               'complex64': PixelsTypecomplex}

def format_s3_uri(uri, endpoint):
    '''
    Combine endpoint and uri
    '''
    parsed_uri = urlsplit(uri)
    url =  "{0.netloc}".format(parsed_uri)
    if endpoint:
        parsed_endpoint = urlsplit(endpoint)
        endpoint = "{0.netloc}".format(parsed_endpoint)
    else:
        endpoint = AWS_DEFAULT_ENDPOINT
    return "{0.scheme}".format(parsed_uri) + "://" + endpoint + "/" + url + "{0.path}".format(parsed_uri)


def load_array(store, path=None):
    arr = zarr.open(store=store, mode="r", path=path)
    return arr


def load_attrs(store, path=None):
    """
    Load the attrs from the root group or path subgroup
    """
    root = zarr.open(store=store, mode="r", path=path)
    attrs = root.attrs.asdict()
    if "ome" in attrs:
        attrs = attrs["ome"]
    return attrs


def masks_from_labels_nd(
        labels_nd, axes="tcz", label_props=None):
    rois = {}

    colors_by_value = {}
    if "colors" in label_props:
        for color in label_props["colors"]:
            pixel_value = color.get("label-value", None)
            rgba = color.get("rgba", None)
            if pixel_value and rgba and len(rgba) == 4:
                colors_by_value[pixel_value] = rgba

    text_by_value = {}
    if "properties" in label_props:
        for props in label_props["properties"]:
            pixel_value = props.get("label-value", None)
            text = props.get("omero:text", None)
            if pixel_value and text:
                text_by_value[pixel_value] = text

    # For each label value, we create an ROI that
    # contains 2D masks for each time point, channel, and z-slice.
    for i in range(1, labels_nd.max() + 1):
        print("Mask value", i)
        if not np.any(labels_nd == i):
            continue

        masks = []
        bin_img = labels_nd == i

        sizes = {dim: labels_nd.shape[axes.index(dim)] for dim in axes}
        size_t = sizes.get("t", 1)
        size_c = sizes.get("c", 1)
        size_z = sizes.get("z", 1)

        for t in range(size_t):
            for c in range(size_c):
                for z in range(size_z):
                    print("t, c, z", t, c, z)

                    indices = []
                    if "t" in axes:
                        indices.append(t)
                    if "c" in axes:
                        indices.append(c)
                    if "z" in axes:
                        indices.append(z)

                    # indices.append(np.s_[::])
                    # indices.append(np.s_[x:x_max:])

                    # slice down to 2D plane
                    plane = bin_img[tuple(indices)]

                    if not np.any(plane):
                        continue

                    # plane = plane.compute()

                    # Find bounding box to minimise size of mask
                    xmask = plane.sum(0).nonzero()[0]
                    ymask = plane.sum(1).nonzero()[0]
                    print("xmask", xmask, "ymask", ymask)
                    # if any(xmask) and any(ymask):
                    x0 = min(xmask)
                    w = max(xmask) - x0 + 1
                    y0 = min(ymask)
                    h = max(ymask) - y0 + 1
                    print("cropping to x, y, w, h", x0, y0, w, h)
                    submask = plane[y0:(y0 + h), x0:(x0 + w)]

                    mask = MaskI()
                    mask.setBytes(np.packbits(np.asarray(submask, dtype=int)))
                    mask.setWidth(rdouble(w))
                    mask.setHeight(rdouble(h))
                    mask.setX(rdouble(x0))
                    mask.setY(rdouble(y0))

                    if i in colors_by_value:
                        ch = ColorHolder.fromRGBA(*colors_by_value[i])
                        mask.setFillColor(rint(ch.getInt()))
                    if "z" in axes:
                        mask.setTheZ(rint(z))
                    if "c" in axes:
                        mask.setTheC(rint(c))
                    if "t" in axes:
                        mask.setTheT(rint(t))
                    if i in text_by_value:
                        mask.setTextValue(rstring(text_by_value[i]))

                    masks.append(mask)

        rois[i] = masks

    return rois


def rois_from_labels_nd(conn, img, labels_nd, axes="tcz", label_props=None):
    # Text is set on Mask shapes, not ROIs
    rois = masks_from_labels_nd(labels_nd, axes, label_props)

    for label, masks in rois.items():
        if len(masks) > 0:
            create_roi(conn, img=img, shapes=masks)


def create_roi(conn, img, shapes, name=None):
    # create an ROI, link it to Image
    roi = RoiI()
    roi.setImage(omero.model.ImageI(img.id, False))
    if name is not None:
        roi.setName(rstring(name))
    for shape in shapes:
        roi.addShape(shape)
    # Save the ROI (saves any linked shapes too)
    print(f"Save ROI for image {img.getName()}")
    return conn.getUpdateService().saveAndReturnObject(roi)


def parse_image_metadata(store, img_attrs, image_path=None):
    """
    Parse the image metadata
    """
    multiscale_attrs = img_attrs['multiscales'][0]
    array_path = multiscale_attrs["datasets"][0]["path"]
    if image_path is not None:
        array_path = image_path.rstrip("/") + "/" + array_path
    # load .zarray from path to know the dimension
    array_data = load_array(store, array_path)
    sizes = {}
    shape = array_data.shape
    axes = multiscale_attrs.get("axes")
    # Need to check the older version
    if axes:
        for axis, size in zip(axes, shape):
            if isinstance(axis, str):
                sizes[axis] = size  # v0.3
            else:
                sizes[axis["name"]] = size

    pixels_type = array_data.dtype.name
    return sizes, pixels_type


def create_labels(conn, store, image, labels_path):

    """
    Create labels for the image
    """
    label_image = load_attrs(store, labels_path)
    
    axes = label_image["multiscales"][0]["axes"]
    axes_names = [axis["name"] for axis in axes]
    label_props = label_image.get("image-label", None)

    ds_path = label_image["multiscales"][0]["datasets"][0]["path"]
    array_path = f"{labels_path}/{ds_path}/"
    labels_nd = load_array(store, array_path)
    labels_data = labels_nd[slice(None)]
    print("labels_nd", labels_nd)
    print("labels_data", labels_data)
    print("axes_names", axes_names)

    # Create ROIs from the labels
    rois_from_labels_nd(conn, image, labels_data, axes_names, label_props)


def create_image(conn, store, image_attrs, object_name, families, models, args, image_path=None):
    '''
    Create an Image/Pixels object
    '''
    query_service = conn.getQueryService()
    pixels_service = conn.getPixelsService()
    sizes, pixels_type = parse_image_metadata(store, image_attrs, image_path)
    size_t = sizes.get("t", 1)
    size_z = sizes.get("z", 1)
    size_x = sizes.get("x", 1)
    size_y = sizes.get("y", 1)
    # if channels is None or len(channels) != size_c:
    channels = list(range(sizes.get("c", 1)))
    omero_pixels_type = query_service.findByQuery("from PixelsType as p where p.value='%s'" % PIXELS_TYPE[pixels_type], None)
    iid = pixels_service.createImage(size_x, size_y, size_z, size_t, channels, omero_pixels_type, object_name, "", conn.SERVICE_OPTS)
    iid = iid.getValue()

    rnd_def = None
    image = conn.getObject("Image", iid)
    omero_attrs = image_attrs.get('omero', None)
    if omero_attrs is not None:
        set_channel_names(conn, iid, omero_attrs)
        # Check rendering settings
        rnd_def = set_rendering_settings(omero_attrs, pixels_type, image.getPixelsId(), families, models)

    img_obj = image._obj
    set_external_info(img_obj, args, image_path)
    
    # check for labels...
    if args.labels:
        labels_path = "labels/"
        if image_path is not None:
            labels_path = image_path.rstrip("/") + "/" + labels_path
        print("checking for labels at", labels_path)
        try:
            labels_attrs = load_attrs(store, labels_path)
            print("labels_attrs", labels_attrs)
            if "labels" in labels_attrs:
                for pth in labels_attrs["labels"]:
                    create_labels(conn, store, image, f"{labels_path}/{pth}/")
        except FileNotFoundError:
            pass
    return img_obj, rnd_def

def hex_to_rgba(hex_color):
    """
    Converts a hex color code to an RGB array.
    """
    if len(hex_color) == 3:
      hex_color = hex_color[0]*2 + hex_color[1]*2 + hex_color[2]*2
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return [r, g, b]


def get_channels(omero_info):
    '''
    Find the name of the channels if specified
    '''
    channel_names = []
    if omero_info is None:
        return channel_names
    for index, entry in enumerate(omero_info.get('channels', [])):
        channel_names.append(entry.get('label', index))
    return channel_names


def set_channel_names(conn, iid, omero_attrs):
    channel_names = get_channels(omero_attrs)
    if len(channel_names) == 0:
        return
    nameDict = dict((i + 1, name) for i, name in enumerate(channel_names))
    conn.setChannelNames("Image", [iid], nameDict)


def set_rendering_settings(omero_info, pixels_type, pixels_id, families, models):
    '''
    Extract the rendering settings and the channels information
    '''
    if omero_info is None:
        return
    rdefs = omero_info.get('rdefs', None)
    if rdefs is None:
        rdefs = dict()
    rnd_def = omero.model.RenderingDefI()
    rnd_def.version = rint(0)
    rnd_def.defaultZ = rint(rdefs.get('defaultZ', 0))
    rnd_def.defaultT = rint(rdefs.get('defaultT', 0))
    value = rdefs.get('model', 'rgb')
    if value == 'color':
        value = 'rgb'
    ref_model = None
    for m in models:
        mv = m.getValue()._val
        if mv == 'rgb':
            ref_model = m
        if mv == value:
            rnd_def.model = m
    if rnd_def.model is None:
        rnd_def.model = ref_model

    q_def = omero.model.QuantumDefI()
    q_def.cdStart = rint(0)
    q_def.cdEnd = rint(255)
    # Flag to select a 8-bit depth (<i>=2^8-1</i>) output interval
    q_def.bitResolution = rint(255)
    rnd_def.quantization = q_def
    rnd_def.pixels = omero.model.PixelsI(pixels_id, False)

    if pixels_type.startswith('float'):
        pixels_min = finfo(pixels_type).min
        pixels_max = finfo(pixels_type).max
    else:
        pixels_min = iinfo(pixels_type).min
        pixels_max = iinfo(pixels_type).max
    for entry in omero_info.get('channels', []):
        cb = omero.model.ChannelBindingI()
        rnd_def.addChannelBinding(cb)
        cb.coefficient = rdouble(entry.get('coefficient', 1.0))
        cb.active = rbool(entry.get('active', False))
        value = entry.get('family', "linear")
        ref_family = None
        for f in families:
            fv = f.getValue()._val
            if fv == "linear":
                ref_family = f
            if fv == value:
                cb.family = f
        if cb.family is None:
            cb.family = ref_family

        # convert color to rgba
        rgb = hex_to_rgba(entry.get('color', "000000").lstrip("#")) # default to black is no color set
        cb.red = rint(rgb[0])
        cb.green = rint(rgb[1])
        cb.blue = rint(rgb[2])
        cb.alpha = rint(255)
        cb.noiseReduction = rbool(False)

        window = entry.get("window", None)
        if window:
            cb.inputStart = rdouble(window.get("start", pixels_min))
            cb.inputEnd = rdouble(window.get("end", pixels_max))
        inverted = entry.get("inverted", False)
        if inverted: # add codomain
            ric = omero.model.ReverseIntensityContextI()
            ric.reverse = rbool(inverted)
            cb.addCodomainMapContext(ric)
    return rnd_def


def load_families(query_service):
    ctx = {'omero.group': '-1'}
    return query_service.findAllByQuery('select f from Family as f', None, ctx)


def load_models(query_service):
    ctx = {'omero.group': '-1'}
    return query_service.findAllByQuery('select f from RenderingModel as f', None, ctx)


def register_image(conn, store, args, img_attrs=None, image_path=None):
    """
    Register the ome.zarr image in OMERO.
    """

    update_service = conn.getUpdateService()
    query_service = conn.getQueryService()
    families = load_families(query_service)
    models = load_models(query_service)

    if img_attrs is None:
        img_attrs = load_attrs(store, image_path)
    if args.name:
        image_name = args.name
    elif "name" in img_attrs:
        image_name = img_attrs["name"]
    else:
        image_name = args.uri.rstrip("/").split("/")[-1]
        if image_path is not None:
            image_name = f"{image_name} [{image_path}]"
    image, rnd_def = create_image(conn, store, img_attrs, image_name, families, models, args, image_path=image_path)
    update_service.saveAndReturnObject(image)
    if rnd_def is not None:
        update_service.saveAndReturnObject(rnd_def)

    print("Created Image", image.id.val)
    return image


def determine_naming(values):
    '''
    Determine the name of columns or rows of a plate
    '''
    if len(values) > 0:
        value = values[0]['name']
        if value.isdigit():
            return "number"
    return "letter"

def create_plate_acquisition(pa):
    '''
    Create a plate acquisition object
    '''
    plate_acquisition = omero.model.PlateAcquisitionI()
    if pa.get("name"):
        plate_acquisition.name = rstring(pa.get("name"))
    else:
        plate_acquisition.name = rstring(pa.get("id"))
    if pa.get("maximumfieldcount"):
        plate_acquisition.maximumFieldCount = rint(pa.get("maximumfieldcount"))
    if pa.get("starttime"):
        plate_acquisition.startTime = rtime(pa.get("starttime"))
    if pa.get("endtime"):
        plate_acquisition.endTime = rtime(pa.get("endtime"))
    return plate_acquisition


# def register_plate(conn, uri, name=None, transport_params=None, endpoint=None, uri_parameters=None):
def register_plate(conn, store, args, attrs):
    '''
    Register a plate
    '''

    plate_attrs = attrs["plate"]

    object_name = args.name
    if object_name is None:
        object_name = plate_attrs.get("name", None)
    if object_name is None:
        object_name = args.uri.rstrip("/").split("/")[-1].split(".")[0]

    update_service = conn.getUpdateService()
    query_service = conn.getQueryService()
    families = load_families(query_service)
    models = load_models(query_service)

    # Create a plate
    plate = omero.model.PlateI()
    plate.name = rstring(object_name)
    plate.columnNamingConvention = rstring(determine_naming(plate_attrs['columns']))
    plate.rowNamingConvention = rstring(determine_naming(plate_attrs['rows']))
    plate.rows = rint(len(plate_attrs['rows']))
    plate.columns = rint(len(plate_attrs['columns']))

    acquisitions = plate_attrs.get('acquisitions')
    plate_acquisitions = {}
    if acquisitions is not None and len(acquisitions) > 1:
        for pa in acquisitions:
            plate_acquisition = create_plate_acquisition(pa)
            plate.addPlateAcquisition(plate_acquisition)

    plate = update_service.saveAndReturnObject(plate)
    print("Plate created with id:", plate.id.val)

    # load the new plate acquisitions and map them to the original IDs
    if acquisitions is not None and len(acquisitions) > 1:
        pwrapper = conn.getObject("Plate", plate.id.val)
        pas = list(pwrapper.listPlateAcquisitions())
        for pa, saved in zip(acquisitions, pas):
            plate_acquisitions[pa["id"]] = saved.id
        print('plate_acquisitions', plate_acquisitions)

    # for bug in omero-cli-zarr - need to handle dupliate Wells!
    well_paths = []

    well_count = len(plate_attrs["wells"])
    for well_index, well_attrs in enumerate(plate_attrs["wells"]):
        images_to_save = []
        rnd_defs = []
        # read metadata
        row_index = well_attrs["rowIndex"]
        column_index = well_attrs["columnIndex"]
        well_path = well_attrs['path']
        if well_path in well_paths:
            continue
        else:
            well_paths.append(well_path)
        print("well_path", well_path, f"({well_index}/{well_count})")
        # create OMERO object
        well = omero.model.WellI()
        well.plate = omero.model.PlateI(plate.getId(), False)
        well.column = rint(column_index)
        well.row = rint(row_index)

        well_attrs = load_attrs(store, well_path)
        well_samples_attrs = well_attrs["well"]["images"]


        for sample_attrs in well_samples_attrs:
            image_path = f"{well_path}/{sample_attrs['path']}/"

            img_attrs = load_attrs(store, image_path)
            image_name = img_attrs.get('name', f"{well_path}/{sample_attrs['path']}")

            image, rnd_def = create_image(conn, store, img_attrs, image_name, families, models, args, image_path)

            images_to_save.append(image)
            if rnd_def is not None:
                rnd_defs.append(rnd_def)
            # Link well sample and plate acquisition
            ws = omero.model.WellSampleI()
            if 'acquisition' in sample_attrs:
                acquisition_id = sample_attrs['acquisition']
                pa_id = plate_acquisitions.get(acquisition_id)
                if pa is not None:
                    ws.plateAcquisition = omero.model.PlateAcquisitionI(pa_id, False)
            ws.image = omero.model.ImageI(image.id.val, False)
            ws.well = well
            well.addWellSample(ws)

        # Save each Well and Images as we go...
        update_service.saveObject(well)
        update_service.saveAndReturnArray(images_to_save)
        if len(rnd_defs) > 0:
            update_service.saveAndReturnIds(rnd_defs)

    print("Plate created with id:", plate.id.val)
    return plate


def set_external_info(image, args, image_path=None):
    '''
    Create the external info and link it to the image
    '''
    extinfo = ExternalInfoI()
    # non-nullable properties
    setattr(extinfo, "entityId", rlong(3))
    setattr(extinfo, "entityType", rstring("com.glencoesoftware.ngff:multiscales"))

    uri = args.uri
    endpoint = args.endpoint
    nosignrequest = args.nosignrequest

    if image_path is not None:
        uri = uri.rstrip("/") + "/" + image_path
    parsed_uri = urlsplit(uri)
    scheme = "{0.scheme}".format(parsed_uri)

    if "http" in scheme:
        endpoint = "https://" + "{0.netloc}".format(parsed_uri)
        nosignrequest = True
        path = "{0.path}".format(parsed_uri)
        if path.startswith("/"):
            path = path[1:]
        uri = "s3://" + path
    
    if not uri.startswith("/"):
        uri = format_s3_uri(uri, endpoint)
    if nosignrequest:
        if not uri.endswith("/"):
            uri = uri + "/"
        uri = uri + "?anonymous=true"
    setattr(extinfo, "lsid", rstring(uri))
    print("lsid:", uri)
    image.details.externalInfo = extinfo

def validate_uri(uri):
    '''
    Check that the protocol is valid and the URI ends with "/"
    '''
    parsed_uri = urlsplit(uri)
    scheme =  "{0.scheme}".format(parsed_uri)
    if "s3" not in scheme:
        raise Exception("Protocol should be s3. Protocol specified is: " + scheme)
    # Check if ends with / otherwise add one
    path = "{0.path}".format(parsed_uri)
    if path.endswith("/"):
        return uri
    return uri + "/"

def validate_endpoint(endpoint):
    '''
    Check that the protocol is valid
    '''
    if endpoint is None or not endpoint:
        return
    parsed_endpoint = urlsplit(endpoint)
    scheme =  "{0.scheme}".format(parsed_endpoint)
    if "https" not in scheme:
        raise Exception("Protocol should be https. Protocol specified is: " + scheme)

def get_uri_parameters(transport_params, nosignrequest):
    if transport_params is None:
        return None
    if nosignrequest:
        return "?anonymous=true"
    return None

def link_to_target(args, conn, obj):
    is_plate = isinstance(obj, omero.model.PlateI)

    if args.target:
        if is_plate:
            target = conn.getObject("Screen", attributes={'id': int(args.target)})
        else:
            target = conn.getObject("Dataset", attributes={'id': int(args.target)})
    else:
        if is_plate:
            target = conn.getObject("Screen", attributes={'name': args.target_by_name})
        else:
            target = conn.getObject("Dataset", attributes={'name': args.target_by_name})

    if target is None:
        print("Target not found")
        return

    if is_plate:
        link = omero.model.ScreenPlateLinkI()
        link.parent = omero.model.ScreenI(target.getId(), False)
        link.child = omero.model.PlateI(obj.getId(), False)
        conn.getUpdateService().saveObject(link)
        print("Linked to Screen", target.getId())
    else:
        link = omero.model.DatasetImageLinkI()
        link.parent = omero.model.DatasetI(target.getId(), False)
        link.child = omero.model.ImageI(obj.getId(), False)
        conn.getUpdateService().saveObject(link)
        print("Linked to Dataset", target.getId())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("uri", type=str, help="The URI to the S3 store")
    parser.add_argument("--endpoint", required=False, type=str, help="Enter the URL endpoint if applicable")
    parser.add_argument("--name", required=False, type=str, help="The name of the image/plate")
    parser.add_argument("--nosignrequest", required=False, action='store_true', help="Indicate to sign anonymously")
    parser.add_argument("--target", required=False, type=str, help="The id of the target (dataset/screen)")
    parser.add_argument("--target-by-name", required=False, type=str, help="The name of the target (dataset/screen)")
    parser.add_argument("--labels", required=False, action='store_true', help="Also import any OME-Zarr labels found")
    
    args = parser.parse_args()

    with cli_login() as cli:
        conn = BlitzGateway(client_obj=cli._client)
        uri = args.uri
        endpoint = args.endpoint
        nosignrequest = args.nosignrequest
        validate_endpoint(endpoint)
        store = None
        if uri.startswith("/"):
            store = zarr.storage.LocalStore(uri, read_only=True)
        else:
            storage_options={}
            if nosignrequest:
                storage_options['anon'] = True

            if endpoint:
                storage_options['client_kwargs'] = {'endpoint_url': endpoint}

            store = zarr.storage.FsspecStore.from_url(uri,
                read_only=True,
                storage_options=storage_options
            )

        zattrs = load_attrs(store)
        objs = []
        if "plate" in zattrs:
            print("Registering: Plate")
            objs = [register_plate(conn, store, args, zattrs)]
        else:
            if "bioformats2raw.layout" in zattrs and zattrs["bioformats2raw.layout"] == 3:
                print("Registering: bioformats2raw.layout")
                series = 0
                series_exists = True
                while series_exists:
                    try:
                        print("Checking for series:", series)
                        obj = register_image(conn, store, args, None, image_path=str(series))
                        objs.append(obj)
                    except FileNotFoundError:
                        series_exists = False
                    series += 1
            else:
                print("Registering: Image")
                objs = [register_image(conn, store, args, zattrs)]

        if args.target or args.target_by_name:
            for obj in objs:
                link_to_target(args, conn, obj)

if __name__ == "__main__":
    main()

