#!/usr/bin/env python
# encoding: utf-8
"""
predictor.py

Created by Elisa Cilia on 2014-08-20.
Copyright (c) 2014 Elisa Cilia. All rights reserved.
"""

import sys
import shutil
import os
import zipfile
from Bio import SeqIO
import re

from . config import *
from . submit import DynaMineJSONInterface
from urllib.request import urlopen
import urllib


class DynaMine:
    def __init__(self, jobdir, light=False, allinoneout=False):
        self._light = light
        self._allinone = allinoneout
        self._jobdir = jobdir

    @staticmethod
    def version():
        version = ''
        try:
            data = urlopen(os.path.join(configuration_file['server_name'], 'version')).read()
            version = 'DynaMine v ' + data.read().strip()
        except Exception as e:
            version = 'ERROR while communicating with the DynaMine server.'
        return version

    def predict(self, filename, output_dir):
        seqs = self._getinputsequences(filename)
        ji = DynaMineJSONInterface(configuration_file['json_api_key'])
        results = ji.submit_sequences(seqs, predictions_only = self._light)
        # print(results)
        if 'status' in results and results['status'] == 'error':
            print(results['message'])
            return False
        else:
            self._save(results, output_dir)
            return True

    def _save(self, results, output_dir):
        filename = results['url'].split('/')[-1]
        path = os.path.join(self._jobdir, filename)
        try:
            with urllib.request.urlopen(results['url']) as dl_file:
                with open(path, 'wb') as out_file:
                    out_file.write(dl_file.read())
            #data = urllib.request.urlopen(results['url'])
            #fd = open(path, 'w')
            #fd.write(data)
            #fd.close()
        except urllib.error.HTTPError:
            sys.stderr.write("\nData %s not found on the server.\n" % (filename))
        except Exception as e:
            sys.stderr.write(str(e))
        with zipfile.ZipFile(path) as zf:
            zf.extractall(self._jobdir)
            output_dir_full_path = f"{self._jobdir}/{output_dir}"
            if (os.path.exists(output_dir_full_path)):
                shutil.rmtree(output_dir_full_path)
            os.rename(f"{self._jobdir}/{filename[:-4]}", output_dir_full_path)
        os.remove(path)

    def _getinputsequences(self, filename):
        seqs = {}
        sequences = SeqIO.parse(filename, "fasta")
        for record in sequences:
            s = re.sub('[^0-9a-zA-Z-]+', '_', record.id.lstrip('>'))
            seqs[s] = str(record.seq)
        return seqs

