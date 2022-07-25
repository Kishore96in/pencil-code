"""
This script assumes pc_start has already been run. Otherwise, we would need to hardcode the values of nx,ny,nz.
"""

import pencil as pc
import numpy as np

dim = pc.read.dim()

a = np.zeros((3,dim.nx,dim.ny,dim.nz))

a[1,:,:,:] = 0.1
a[2,:,:,:] = 0.2
a[0,1,1,1] = 0.05
a[0,2,1,0] = 0.01

pc.util.write_forcing_cont(a)
