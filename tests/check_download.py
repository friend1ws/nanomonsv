#! /usr/bin/env python

import os
import requests

def check_download(url, output_path):

    if not os.path.exists(os.path.dirname(output_path)):
       os.makedirs(os.path.dirname(output_path))

    if not os.path.exists(output_path):
        response = requests.get(url)
        with open(output_path, 'wb') as hout:
            hout.write(response.content)

