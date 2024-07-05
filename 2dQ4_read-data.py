#!/usr/bin/python3
import os
import sys
import h5py
# ------------------------------------------------------------------
# Cycle over 2d Q=4 SYM ensembles in HDF5 file,
# printing group path and quantities of interest

# Annoyingly, need to use absolute Barkla path
path = "/mnt/lustre/users/anosh/susy_bound/susy_bound/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)

def visit_func(name, node):
  if isinstance(node, h5py.Group) and 'zeta' in name:
    print(name, end=' ')
#    print(node.attrs['Ntraj'])
#    print(node.attrs['acceptance'])
#    print(node.attrs['autocorrelation time (|ML|)'], \
#          node.attrs['autocorrelation time (eig)'], \
#          node.attrs['autocorrelation time (bilin)'], \
#          node.attrs['autocorrelation time (Tr[X^2])'], \
#          node.attrs['block size'], node.attrs['Nblocks'])
#    dset = node.get('Maldacena_loop')
#    print(dset.attrs['magnitude ave'], dset.attrs['magnitude err'])
#    dset = node.get('exp_dS')
#    print(dset.attrs['ave'], dset.attrs['err'])
    dset = node.get('SB')
    print(dset.attrs['ave'], dset.attrs['err'], end=' ')
    dset = node.get('bilin')
    print(dset.attrs['ave'], dset.attrs['err'])

with h5py.File("2dQ4_data.h5", 'r') as f:  # Don't need f.close()
  f.visititems(visit_func)
# ------------------------------------------------------------------
