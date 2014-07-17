#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print single Wilson loop measurement as function of 3d-projected A4* r
# rather than (x, y, z) displacements

# Parse argument: the file to consider
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])

# Convenience constants
TOL = 1.0e-6
invSq2  = 1.0 / np.sqrt(2)
invSq6  = 1.0 / np.sqrt(6)
invSq12 = 1.0 / np.sqrt(12)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Need to construct lists of scalar distances and their multiplicities
# while we accumulate the corresponding data
# This list has three components: r, its multiplicity and the running sum
r = []
for line in open(toOpen):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('hvy_pot: MAX_T '):
    MAX_T = int((line.split())[3].rstrip(','))    # Strip ',' from end
  # Format: POT_LOOP x y z t dat      (similarly for D_LOOP and POLAR_LOOP)
  elif line.startswith('POT_LOOP '):
    temp = line.split()
    if int(temp[4]) > 1:
      break                 # Only consider a single value of t=1!
    x = float(temp[1])
    y = float(temp[2])
    z = float(temp[3])
    dat = float(temp[5])
    x_a4 = (x - y) * invSq2
    y_a4 = (x + y - 2.0 * z) * invSq6
    z_a4 = (x + y + z) * invSq12
    this_r = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2)

    # Check all possible periodic shifts to find true r
    for xx in [x + L, x, x - L]:
      for yy in [y + L, y, y - L]:
        for zz in [z + L, z, z - L]:
          xx_a4 = (xx - yy) * invSq2
          yy_a4 = (xx + yy - 2.0 * zz) * invSq6
          zz_a4 = (xx + yy + zz) * invSq12
          test_r = np.sqrt(xx_a4**2 + yy_a4**2 + zz_a4**2)
          if test_r - this_r < -TOL:   # Avoid negative roundoff
#            print "|(%d, %d, %d)| = %.6g -->" \
#                  % (int(x), int(y), int(z), this_r),
#            print "|(%d, %d, %d)| = %.6g" \
#                  % (int(xx), int(yy), int(zz), test_r)
            this_r = test_r
    print x, y, z, "-->", this_r

    # Inefficient, but whatever
    done = -1
    for i in range(len(r)):
      if np.fabs(r[i][0] - this_r) < TOL:
        r[i][1] += 1
        r[i][2] += dat
        done = 1
        break
    if done < 0:
      r.append([this_r, 1, dat])
Npts = len(r)

# Sort by magnitude (column zero), not count
r = sorted(r, key=lambda x: x[0])
for i in range(Npts):   # Print r and average
  print "%.4g %.4g" % (r[i][0], r[i][2]  / r[i][1])
# ------------------------------------------------------------------

