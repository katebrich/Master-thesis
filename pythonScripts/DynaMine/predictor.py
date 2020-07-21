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

    def predict(self, fasta_path, pdb_id, chain_id):
        seqs = self._getinputsequences(fasta_path, pdb_id, chain_id)
        ji = DynaMineJSONInterface(configuration_file['json_api_key'])
        results = ji.submit_sequences(seqs, predictions_only = self._light)
        # print(results)
        if 'status' in results and results['status'] == 'error':
            print(results['message'])
            return False
        else:
            self._save(results, pdb_id, chain_id)
            return True

    def _save(self, results, pdb_id, chain_id):
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
            output_full_path = f"{self._jobdir}/{pdb_id}{chain_id}.txt"
            unzipped_path = f"{self._jobdir}/{filename[:-4]}"
            if (os.path.exists(output_full_path)):
                os.remove(output_full_path)
            os.rename(f"{unzipped_path}/{pdb_id}{chain_id}_backbone.pred", output_full_path)
        os.remove(path)
        shutil.rmtree(unzipped_path)

    def _getinputsequences(self, fasta_path, pdb_id, chain_id):
        seqs = {}
        with open(fasta_path, 'r') as f:
            seq = f.read()
            seqs[f"{pdb_id}{chain_id}"] = seq
        return seqs

