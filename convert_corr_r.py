#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print single SUGRA measurement as function of A4* scalar distances r
# rather than (x, y, z, t) displacements
# Ignore Konishi for now

# Parse argument: the file to consider
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])

# Convenience constants for sanity check
TOL = 1.0e-6
invSq2  = 1.0 / np.sqrt(2)
invSq6  = 1.0 / np.sqrt(6)
invSq12 = 1.0 / np.sqrt(12)
invSq20 = 1.0 / np.sqrt(20)
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
  if line.startswith('d_correlator_r: MAX_T = '):
    MAX_X = int((line.split())[6])
  # Format: CORR_S x y z t dat
  elif line.startswith('CORR_S '):
    temp = line.split()
    x = float(temp[1])
    y = float(temp[2])
    z = float(temp[3])
    t = float(temp[4])
    dat = float(temp[5])
    if t > MAX_X:         # Drop large separations
      continue
    this_r = np.sqrt((x**2 + y**2 + z**2 + t**2) * 0.8 \
                     - 2.0 * (x * (y + z + t) + y * (z + t) + z * t) * 0.2)
#    print x, y, z, t, "-->", this_r

    # Sanity check
    x_a4 = (x - y) * invSq2
    y_a4 = (x + y - 2.0 * z) * invSq6
    z_a4 = (x + y + z - 3.0 * t) * invSq12
    t_a4 = (x + y + z + t) * invSq20
    check_r = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2 + t_a4**2)
    if np.fabs(this_r - check_r) > TOL:
      print "ERROR:", this_r, "isn't", check_r
      sys.exit(1)

    # Wrapping check


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

