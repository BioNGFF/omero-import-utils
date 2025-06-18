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
import os

import argparse
import boto3
import botocore
import botocore.client
import json
from smart_open import open
from numpy import dtype, iinfo

from getpass import getpass

from omero.gateway import BlitzGateway
from omero.gateway import OMERO_NUMPY_TYPES

import omero

from omero.model.enums import PixelsTypeint8, PixelsTypeuint8, PixelsTypeint16
from omero.model.enums import PixelsTypeuint16, PixelsTypeint32
from omero.model.enums import PixelsTypeuint32, PixelsTypefloat
from omero.model.enums import PixelsTypecomplex, PixelsTypedouble

from omero.model import ExternalInfoI
from omero.rtypes import rbool, rdouble, rint, rlong, rstring


EXTENSION_JSON = "zarr.json"

ENDPOINT = "s3.amazonaws.com"

OBJECT_PLATE = "plate"
OBJECT_IMAGE = "image"

PIXELS_TYPE = {'int8': PixelsTypeint8,
               'int16': PixelsTypeint16,
               'uint8': PixelsTypeuint8,
               'uint16': PixelsTypeuint16,
               'int32': PixelsTypeint32,
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
        endpoint = ENDPOINT
    return "{0.scheme}".format(parsed_uri) + "://" + endpoint + "/" + url + "{0.path}".format(parsed_uri)


def create_client(endpoint, nosignrequest=False):
    """
    Create a boto3 client to connect to S3
    """
    config = None
    if nosignrequest:
        config = botocore.client.Config(signature_version=botocore.UNSIGNED)
    session = boto3.Session()
    if endpoint:
        if config:
            client = session.client('s3', endpoint_url=endpoint, config=config)
        else:
            client = session.client('s3', endpoint_url=endpoint)
    else:
        if config:
            client = session.client('s3', config=config)
        else:
            client = session.client('s3')
    transport_params = {'client': client}
    return transport_params

def determine_object_to_register(uri, transport_params):
    """
    Determine the object to register: supported Plate and Image
    """
    extension = ".zattrs"
    path = uri + extension
    try:
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)
    except Exception as e:
        extension = EXTENSION_JSON
        path = uri + extension
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)

    if extension == EXTENSION_JSON:
        if "plate" in zattrs['attributes']['ome']:
            return OBJECT_PLATE, uri
    else:
        if "plate" in zattrs:
            return OBJECT_PLATE, uri
    if "bioformats2raw.layout" in zattrs and zattrs["bioformats2raw.layout"] == 3:
        uri = f"{uri}0/"
    return OBJECT_IMAGE, uri


def parse_image_metadata(uri, image_zattrs, transport_params, extension):
    """
    Parse the image metadata
    """
    array_path = image_zattrs["datasets"][0]["path"]
    # load .zarray from path to know the dimension
    if extension == EXTENSION_JSON:
        data_type_key = "data_type"
        path = uri + array_path + "/zarr.json"
    else:
        data_type_key = "dtype"
        path = uri + array_path + "/.zarray"

    with open(path, 'rb', transport_params=transport_params) as f:
        array_data = json.load(f)
    sizes = {}
    shape = array_data["shape"]
    axes = image_zattrs.get("axes")
    # Need to check the older version
    if axes:
        for dim, size in zip(axes, shape):
            sizes[dim["name"]] = size

    pixels_type = dtype(array_data[data_type_key]).name
    return sizes, pixels_type


def create_image(pixels_service, query_service, sizes, pixels_type, object_name, options):
    '''
    Create an Image/Pixels object
    '''
    size_t = sizes.get("t", 1)
    size_z = sizes.get("z", 1)
    size_x = sizes.get("x", 1)
    size_y = sizes.get("y", 1)
    channels = list(range(sizes.get("c", 1)))
    omero_pixels_type = query_service.findByQuery("from PixelsType as p where p.value='%s'" % PIXELS_TYPE[pixels_type], None)
    return pixels_service.createImage(size_x, size_y, size_z, size_t, channels, omero_pixels_type,  object_name, "", options)

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


def set_rendering_settings(omero_info, pixels_type, pixels_id, families, models):
    '''
    Extract the rendering settings and the channels information
    '''
    if omero_info is None:
        return
    channel_names = {}
    rdefs = omero_info.get('rdefs', None)
    rnd_def = None
    if rdefs is not None:
        rnd_def = omero.model.RenderingDefI()
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

    pixels_min = iinfo(pixels_type).min
    pixels_max = iinfo(pixels_type).max
    for index, entry in enumerate(omero_info.get('channels', [])):
        channel_names[index] = entry.get('label', index)
        if rnd_def is not None:
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
            rgb = hex_to_rgba(entry.get('color', "000000")) # default to black is no color set
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
    params = omero.sys.ParametersI()
    models = query_service.findAllByQuery('select f from Family as f', params, ctx)


def load_models(query_service):
    ctx = {'omero.group': '-1'}
    params = omero.sys.ParametersI()
    models = query_service.findAllByQuery('select f from RenderingModel as f', params, ctx)
 

def register_image(uri, transport_params,  host, username, password, name="", endpoint=""):
    """
    Register the ome.zarr image in OMERO.
    """
    extension = ".zattrs"
    try:
        path = uri + extension
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)
        image_zattrs = zattrs['multiscales'][0]
        omero_zattrs = zattrs.get('omero', None)
    except Exception as e:
        extension = "zarr.json"
        path = uri + extension
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)
        omero_zattrs = zattrs["attributes"]["ome"].get('omero', None)
        image_zattrs = zattrs["attributes"]["ome"]['multiscales'][0]

    sizes, pixels_type = parse_image_metadata(uri, image_zattrs, transport_params, extension)
    if name:
        object_name = name
    else:
        object_name = uri.rstrip("/").split("/")[-1]
    # connect to omero
    try:
        conn = BlitzGateway(username, password, host=host, secure=True)
        conn.connect()
        conn.c.enableKeepAlive(60)
        session = conn.c.getSession()
        pixels_service = session.getPixelsService()
        query_service = session.getQueryService()
        update_service = session.getUpdateService()
        iid = create_image(pixels_service, query_service, sizes, pixels_type, object_name, conn.SERVICE_OPTS)
        # load the image object
        image = conn.getObject("Image", iid)
        # register external info
        set_external_info(uri, endpoint, image)
        # save the image with the updated external infor
        image.save()
        # Set the rendering object
        # loads the families and the color model
        ctx = {'omero.group': '-1'}
        params = omero.sys.ParametersI()
        families = query_service.findAllByQuery('select f from Family as f', params, ctx)
        models = query_service.findAllByQuery('select f from RenderingModel as f', params, ctx)
        rnd_def = set_rendering_settings(omero_zattrs, pixels_type, image.getPixelsId(), families, models)
        # Save the rendering def
        update_service.saveAndReturnObject(rnd_def)

    except Exception as e:
        raise e
    finally:
        conn.close()

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
        plate_acquisition.startTime = rint(pa.get("starttime"))
    if pa.get("endtime"):
        plate_acquisition.endTime = rint(pa.get("endtime"))
    return plate_acquisition
    

def register_plate(uri, transport_params, host, username, password, name="", endpoint=""):
    '''
    Register a plate
    '''
    extension = ".zattrs"
    path = uri + extension
    try:
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)
    except Exception as e:
        extension = "zarr.json"
        path = uri + extension
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)
    

    if extension == EXTENSION_JSON:
        plate_attrs = zattrs["attributes"]["ome"]["plate"]
    else:
        plate_attrs = zattrs["plate"]

    if name:
        object_name = name
    else:
        object_name = uri.rstrip("/").split("/")[-1].split(".")[0]

    # connect to omero

    try:
        conn = BlitzGateway(username, password, host=host, secure=True)
        conn.connect()
        conn.c.enableKeepAlive(60)
        session = conn.c.getSession()
        pixels_service = session.getPixelsService()
        query_service = session.getQueryService()
        update_service = session.getUpdateService()

        # loads the families and the color model
        ctx = {'omero.group': '-1'}
        params = omero.sys.ParametersI()
        families = query_service.findAllByQuery('select f from Family as f', params, ctx)
        models = query_service.findAllByQuery('select f from RenderingModel as f', params, ctx)
       
        # Create a plate
        plate = omero.model.PlateI()
        plate.name = rstring(object_name)
        plate.columnNamingConvention = rstring(determine_naming(plate_attrs['columns']))
        plate.rowNamingConvention = rstring(determine_naming(plate_attrs['rows']))
        plate.rows = rint(len(plate_attrs['rows']))
        plate.columns = rint(len(plate_attrs['columns']))
        
        acquisitions = plate_attrs['acquisitions']
        plate_acquisitions = {}
        if len(acquisitions) > 1:
            for pa in acquisitions:
                plate_acquisition =  update_service.saveAndReturnObject(create_plate_acquisition(pa))
                plate_acquisitions[pa.get("id")] = plate_acquisition
                plate.addPlateAcquisition(omero.model.PlateAcquisitionI(plate_acquisition.getId(), False))

        plate = update_service.saveAndReturnObject(plate)
        wells = []
        rnd_defs = []
        images_to_save = []
        for well_attrs in plate_attrs["wells"]:
        	# read metadata
            row_index = well_attrs["rowIndex"]
            column_index = well_attrs["columnIndex"]
            well_path = well_attrs['path']
            # create OMERO object
            well = omero.model.WellI()
            well.plate = omero.model.PlateI(plate.getId(), False)
            well.column = rint(column_index)
            well.row = rint(row_index)

            well_full_path = f"{uri}{well_path}/{extension}"
            with open(well_full_path, 'rb', transport_params=transport_params) as wfp:
                well_samples_json = json.load(wfp)
            if extension == EXTENSION_JSON:
                well_samples_attrs = well_samples_json["attributes"]["ome"]["well"]["images"]
            else:
                well_samples_attrs = well_samples_json["well"]["images"]

            for index, sample_attrs in enumerate(well_samples_attrs):
                image_uri = f"{uri}{well_path}/{sample_attrs['path']}/"

                with open(image_uri + EXTENSION_JSON, 'rb', transport_params=transport_params) as ifp:
                    image_json = json.load(ifp)

                if extension == EXTENSION_JSON:
                    image_attrs = image_json["attributes"]["ome"]['multiscales'][0]
                    omero_attrs = image_json["attributes"]["ome"].get("omero", None)
                else:
                    image_attrs = image_json['multiscales'][0]
                    omero_attrs = image_json.get("omero", None)

                sizes, pixels_type = parse_image_metadata(image_uri, image_attrs, transport_params, extension)
                image_path = f"{image_uri}{image_attrs['datasets'][0]['path']}/"
                iid = create_image(pixels_service, query_service, sizes, pixels_type, image_attrs['name'], conn.SERVICE_OPTS)

                
                # load the image object
                image = conn.getObject("Image", iid)
                set_external_info(image_path, endpoint, image)
                images_to_save.append(image)
                # Check  rendering settings
                rnd_def = set_rendering_settings(omero_attrs, pixels_type, image.getPixelsId(), families, models)
                rnd_defs.append(rnd_def)
                # Link well sample and plate acquisition
                ws = omero.model.WellSampleI()
                acquisition_id = sample_attrs['acquisition']
                pa = plate_acquisitions.get(acquisition_id)
                if pa is not None:
                    ws.plateAcquisition = omero.model.PlateAcquisitionI(pa.getId(), False)
                ws.image = omero.model.ImageI(iid, False)
                ws.well = well
                well.addWellSample(ws)
            #update_service.saveObject(well)

            wells.append(well)
        # Save the 3 arrays
        update_service.saveAndReturnIds(wells)
        update_service.saveAndReturnIds(rnd_defs)
        update_service.saveAndReturnIds(images_to_save)
    except Exception as e:
        raise e
    finally:
        conn.close()

def set_external_info(uri, endpoint, image):
    '''
    Create the external info and link it to the image
    '''
    extinfo = ExternalInfoI()
    # non-nullable properties
    setattr(extinfo, "entityId", rlong(3))
    setattr(extinfo, "entityType", rstring("com.glencoesoftware.ngff:multiscales"))
    setattr(extinfo, "lsid", rstring(format_s3_uri(uri, endpoint)))
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("uri", type=str, help="The URI to the S3 store")
    parser.add_argument("--endpoint", required=False, type=str, help="Enter the URL endpoint if applicable")
    parser.add_argument("--name", required=False, type=str, help="The name of the plate")
    parser.add_argument("--nosignrequest", required=False, action='store_true', help="Indicate to sign anonymously")
    host = input("Host [localhost]: ") or 'localhost'  # noqa
    username = input("Username [root]: ") or 'root'
    password = getpass("Password: ")
    args = parser.parse_args()
    uri = args.uri
    endpoint = args.endpoint
    uri = validate_uri(uri)
    validate_endpoint(endpoint)
    transport_params = create_client(endpoint, args.nosignrequest)
    type_to_register, uri = determine_object_to_register(uri, transport_params)
    if type_to_register == OBJECT_PLATE:
    	register_plate(uri, transport_params, host, username, password, name=args.name, endpoint=endpoint)
    else:
        register_image(uri, transport_params, host, username, password, name=args.name, endpoint=endpoint)

if __name__ == "__main__":
    main()
