#!/usr/bin/env python
# encoding: utf-8
"""
config.py

Created by Elisa Cilia on 2014-08-20.
Copyright (c) 2014 Elisa Cilia. All rights reserved.
"""

import os
import string

base_path = os.path.abspath(os.path.dirname(__file__))

configuration_file = {}
f = open(os.path.join(base_path, "./config.txt"), "r")
lines = [line.strip() for line in f.readlines()]
lines = filter(lambda x: x != "", lines)
f.close()
for line in lines:
    pos = str.find(line, "#")
    if pos != -1:
        line = line[:pos]
    if line:
        variable = line.split("=")
        configuration_file[variable[0].strip()] = eval(variable[1])
