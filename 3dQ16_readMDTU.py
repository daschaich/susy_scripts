#!/usr/bin/python3
import os
import sys
import h5py
# ------------------------------------------------------------------
# Cycle over 3d SYM ensembles in HDF5 file, printing group path and MDTU

# Annoyingly, need to use absolute Barkla path
path = "/users/schaich/SYM/3d/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)

def visit_func(name, node):
  if isinstance(node, h5py.Group) and 'zeta' in name:
    print(name, end=' ')
    print(node.attrs['Ntraj'])

with h5py.File("3dSYM_data.h5", 'r') as f:  # Don't need f.close()
  f.visititems(visit_func)
# ------------------------------------------------------------------
