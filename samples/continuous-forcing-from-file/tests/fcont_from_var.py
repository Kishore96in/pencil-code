#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')
import pencil as pc
import numpy as np

"""
Check if the file created by generate_forcing_cont.py was correctly read by forcing.f90
"""

var = pc.read.var(trimall=True, datadir="../data")
force_from_var = np.stack([var.fx, var.fy, var.fz])
with open('fcont_from_var.out', 'w') as f:
	f.write(np.array2string(force_from_var))
