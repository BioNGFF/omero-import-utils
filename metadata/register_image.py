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
from numpy import dtype

from getpass import getpass

from omero.gateway import BlitzGateway
from omero.gateway import OMERO_NUMPY_TYPES

from omero.model.enums import PixelsTypeint8, PixelsTypeuint8, PixelsTypeint16
from omero.model.enums import PixelsTypeuint16, PixelsTypeint32
from omero.model.enums import PixelsTypeuint32, PixelsTypefloat
from omero.model.enums import PixelsTypecomplex, PixelsTypedouble

from omero.model import ExternalInfoI
from omero.rtypes import rstring, rlong


ENDPOINT = "s3.amazonaws.com"

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

def get_s3_uri(uri, endpoint):
    parsed_uri = urlsplit(uri)
    url =  "{0.netloc}".format(parsed_uri)
    if endpoint:
        parsed_endpoint = urlsplit(endpoint)
        endpoint = "{0.netloc}".format(parsed_endpoint)
    else:
        endpoint = ENDPOINT
    return "{0.scheme}".format(parsed_uri) + "://" + endpoint + "/" + url + "{0.path}".format(parsed_uri)

    
def parse_metadata(uri, endpoint, nosignrequest=False):
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
        path = uri + ".zattrs"
        with open(path, 'rb', transport_params=transport_params) as f:
            zattrs = json.load(f)

        array_path = zattrs["multiscales"][0]["datasets"][0]["path"]
        # load .zarray from path to know the dimension
        path = uri + array_path + "/.zarray"

        with open(path, 'rb', transport_params=transport_params) as f:
             array_data = json.load(f)

        object_name = uri.rstrip("/").split("/")[-1]

        sizes = {}
        shape = array_data["shape"]
        axes = zattrs["multiscales"][0].get("axes")
        # Need to check the older version
        if axes is not None:
            for dim, size in zip(axes, shape):
                sizes[dim["name"]] = size

        pixels_type = dtype(array_data["dtype"]).name
        return sizes, pixels_type, object_name


def register_image(uri, endpoint, nosignrequest=False, name="", host="localhost", username="root", password=""):
    sizes, pixels_type, object_name = parse_metadata(uri, endpoint, nosignrequest)
    if name:
        object_name = name
    # connect to omero
    try:
        conn = BlitzGateway(username, password, host=host, secure=True)
        conn.c.enableKeepAlive(60)
        session = conn.c.getSession()
        pixels_service = session.getPixelsService()
        query_service = session.getQueryService()
        size_t = sizes.get("t", 1)
        size_z = sizes.get("z", 1)
        size_x = sizes.get("x", 1)
        size_y = sizes.get("y", 1)
        channels = list(range(sizes.get("c", 1)))
        omero_pixels_type = query_service.findByQuery("from PixelsType as p where p.value='%s'" % PIXELS_TYPE[pixels_type], None)
        iid = pixels_service.createImage(size_x, size_y, size_z, size_t, channels, omero_pixels_type,  object_name, "", conn.SERVICE_OPTS)

        # register external info
        extinfo = ExternalInfoI()
        # non-nullable properties
        setattr(extinfo, "entityId", rlong(3))
        setattr(extinfo, "entityType", rstring("com.glencoesoftware.ngff:multiscales"))
        setattr(extinfo, "lsid", rstring(get_s3_uri(uri, endpoint)))
        image = conn.getObject("Image", iid.getValue())
        image.details.externalInfo = extinfo
        image.save()

    except Exception as e:
        raise e
    finally:
        conn.close()
    



def validate_uri(uri):
    parsed_uri = urlsplit(uri)
    scheme =  "{0.scheme}".format(parsed_uri)
    if "s3" not in scheme:
        raise Exception("Protocol should be s3. Protocol specified is: " + scheme)

def validate_endpoint(endpoint):
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
    parser.add_argument("--name", required=False, type=str, help="The name of the image")
    parser.add_argument("--nosignrequest", required=False, action='store_true', help="Indicate to sign anonymously")
    host = input("Host [localhost]: ") or 'localhost'  # noqa
    username = input("Username [root]: ") or 'root'
    password = getpass("Password: ")
    args = parser.parse_args()
    uri = args.uri
    endpoint = args.endpoint
    validate_uri(uri)
    validate_endpoint(endpoint)
    register_image(uri, endpoint, args.nosignrequest, args.name, host=host, username=username, password=password)

if __name__ == "__main__":
    main()
