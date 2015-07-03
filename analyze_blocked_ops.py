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
# 8nt8 --> 4nt4
#for la in ['1.0']:
#  smallFile = 'Nc2_4nt4/l' + la + '_b0.4_G0.05/results/ops.dat'
#  largeFile = 'Nc2_8nt8/l' + la + '_b0.4_G0.05/results/ops.dat'
# 12nt12 --> 6nt6
for tag in ['l1.0_b0.6_G0.05']:
  smallFile = 'Nc2_6nt6/' + tag + '/results/ops.dat'
  largeFile = 'Nc2_12nt12/' + tag + '/results/ops.dat'

  # First make sure we're calling this from the right place
  # Should be able to retain these independent of commenting out above
  if not os.path.isfile(smallFile):
    print "ERROR: missing file", smallFile
    sys.exit(1)
  if not os.path.isfile(largeFile):
    print "ERROR: missing file", largeFile
    sys.exit(1)

  # Print out rescaling parameter xi for each blocking level
  # Require that n-times-blocked L Konishi operator
  # matches (n-1)-times-blocked (L/2) result
  smallOp  = []
  smallErr = []
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    smallOp.append(float(temp[1]))
    smallErr.append(float(temp[2]) / float(temp[1]))
  largeOp  = []
  largeErr = []
  for line in open(largeFile):
    if line.startswith('# ') or line.startswith('0 '):
      continue
    temp = line.split()
    largeOp.append(float(temp[1]))
    largeErr.append(float(temp[2]) / float(temp[1]))

  blmax = len(smallOp)
  if not len(largeOp) == blmax:
    print "ERROR: picked up different numbers of blockings...",
    print len(largeOp), "vs.", blmax
    sys.exit(1)

  print "0 1.0     0.0"
  for bl in range(blmax):
    xi = smallOp[bl] / largeOp[bl]                    # Presumably xi^4
    xi_err = xi * np.sqrt(smallErr[bl]**2 + largeErr[bl]**2)
    print "%d %.6g %.4g" % (bl + 1, xi, xi_err)
# ------------------------------------------------------------------
