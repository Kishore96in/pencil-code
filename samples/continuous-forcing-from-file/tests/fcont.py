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

#Check if the file created by generate_forcing_cont.py was correctly read by forcing.f90
var = pc.read.var(trimall=True, datadir="../data")
force_from_var = np.stack([var.fx, var.fy, var.fz])
with open('fcont_from_var.out', 'w') as f: #TODO: Does this need to have the same name as the .py script? Presumably it only needs to have extension '.out'.
	f.write(np.array2string(force_from_var))

#Write out a continuous forcing file
b = np.arange(0,81).reshape((3,3,3,3))
pc.util.write_forcing_cont(b, outfile="fcont_write_test.out")
