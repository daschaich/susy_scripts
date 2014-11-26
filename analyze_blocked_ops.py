#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compare 6nt6 and 12nt12 Konishi operators

# No arguments yet...
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
for la in ['0.5', '1.0', '1.5']:
  smallFile = 'Nc2_6nt6/l' + la + '_b0.5_f0.0_k0.5/results/ops.dat'
  largeFile = 'Nc2_12nt12/l' + la + '_b0.5_f0.0_k0.5/results/ops.dat'
  if not os.path.isfile(smallFile) or not os.path.isfile(largeFile):
    print "ERROR: missing file for lambda =", la
    sys.exit(1)

  # Print out rescaling parameter xi for each blocking level
  # Require that n-times-blocked 12nt12 Konishi operator
  # matches (n-1)-times-blocked 6nt6 result
  smallOp  = []
  smallErr = []
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    else:
      temp = line.split()
      smallOp.append(float(temp[1]))
      smallErr.append(float(temp[2]) / float(temp[1]))
  largeOp  = []
  largeErr = []
  for line in open(largeFile):
    if line.startswith('# ') or line.startswith('0 '):
      continue
    else:
      temp = line.split()
      largeOp.append(float(temp[1]))
      largeErr.append(float(temp[2]) / float(temp[1]))

  bl = len(smallOp)
  if not len(largeOp) == bl:
    print "ERROR: picked up different numbers of blockings..."
    sys.exit(1)

  for i in range(bl):
    xi = smallOp[i] / largeOp[i]                    # Presumably xi^4
    xi_err = xi * np.sqrt(smallErr[i]**2 + largeErr[i]**2)
    print "%.4g %d %.6g %.4g" % (float(la), i, xi, xi_err)
# ------------------------------------------------------------------
