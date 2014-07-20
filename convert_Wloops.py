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
# Map (x, y, z) to r on L^3 reduction of A4* lattice
# Check all possible periodic shifts to find true r
def A4map(x_in, y_in, z_in, L):
  r = np.sqrt((x_in**2 + y_in**2 + z_in**2) * 0.75 \
              - (x_in * (y_in + z_in) + y_in * z_in) * 0.5)
  for x in [x_in + L, x_in, x_in - L]:
    for y in [y_in + L, y_in, y_in - L]:
      for z in [z_in + L, z_in, z_in - L]:
        test = np.sqrt((x**2 + y**2 + z**2) * 0.75 \
                       - (x * (y + z) + y * z) * 0.5)

        # Sanity check -- can be commented out for more speed
        x_a4 = (x - y) * invSq2
        y_a4 = (x + y - 2.0 * z) * invSq6
        z_a4 = (x + y + z) * invSq12
        check = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2)
        if np.fabs(test - check) > TOL:
          print "ERROR: %.6g isn't %.6g for (%d, %d, %d)" \
                % (check, test, x, y, z)
          sys.exit(1)

        # True r is the smallest
        # Try to avoid negative roundoff
        # so that we can see when periodic shifts are really necessary
        if test - r < -TOL:
          print "|(%d, %d, %d)| = %.6g"     % (x_in, y_in, z_in, r),
          print "--> |(%d, %d, %d)| = %.6g" % (x, y, z, test)
          r = test
  print "%d %d %d --> %.6g" % (x_in, y_in, z_in, r)
  return r
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
  # Format: hvy_pot: MAX_T = #, MAX_X = #
  elif line.startswith('hvy_pot: MAX_T '):
    temp = line.split()
    MAX_T = int((temp[3]).rstrip(','))    # Strip ',' from end
    MAX_X = int(temp[6])

    # Find smallest r that is cut due to MAX_X
    MAX_r = 100.0 * MAX_X   # To be overwritten
    for y in range(MAX_X + 1):
      for z in range(MAX_X + 1):
        test = A4map(MAX_X + 1, y, z, L)
        if test < MAX_r:
          MAX_r = test
    print "MAX_r = %.6g" % MAX_r
    print "--------------------------------------"

  # Format: POT_LOOP x y z t dat      (similarly for D_LOOP and POLAR_LOOP)
  elif line.startswith('POT_LOOP '):
    temp = line.split()
    if int(temp[4]) > 1:
      break                 # Only consider a single value of t=1!
    x = int(temp[1])
    y = int(temp[2])
    z = int(temp[3])
    dat = float(temp[5])
    this_r = A4map(x, y, z, L)

    # Accumulate multiplicities if 0.5 < this_r < MAX_r
    # Try to avoid roundoff issues in latter comparison
    if this_r < 0.5 or this_r - MAX_r > -TOL:
      continue
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
print "--------------------------------------"
for i in range(Npts):   # Print r and average
  print "%.6g %.6g" % (r[i][0], r[i][2]  / r[i][1])
# ------------------------------------------------------------------

