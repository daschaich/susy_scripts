#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Check that matrix elements <psi_j | D | psi_i> of DDdag eigenvectors
# come in positive/negative pairs
TOL = 1e-6

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Cycle over eigenvalue files
for filename in glob.glob('Out/eig.*'):
  fail = -1
  r = []          # List of real and imaginary parts of the matrix elements
  for line in open(filename):
    # Format: D_eig 0 (re, im)
    if line.startswith('D_eig '):
      temp = line.split('(')
      re = float(((temp[1]).split(','))[0])
      temp = line.split(', ')
      im = float(((temp[1]).split(')'))[0])
      r.append([re, im])

  # Sort by real part and check whether r[i] = -r[N - 1 - i]
  r = sorted(r, key=lambda x: x[0])
  N = len(r)
  for i in range(N / 2):
    for j in range(2):
      if np.fabs(r[i][j] + r[N - 1 - i][j]) > TOL:
        fail = 1
  if fail > 0:
    print "Apparent pairing breakdown for", filename
# ------------------------------------------------------------------
