#! /usr/bin/env python

import os

def check_download(url, output_path):
    
    import requests

    if not os.path.exists(os.path.dirname(output_path)):
       os.makedirs(os.path.dirname(output_path))

    if not os.path.exists(output_path):
        response = requests.get(url)
        with open(output_path, 'wb') as hout:
            hout.write(response.content)

def check_download_s3(bucket_name, object_name, file_name):

    import boto3
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name))

    if not os.path.exists(file_name):
        s3 = boto3.client('s3')
        s3.download_file(bucket_name, object_name, file_name)

