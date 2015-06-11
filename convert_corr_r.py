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
# Map (x, y, z, t) to r on L^3 x Nt A4* lattice
# Check all possible periodic shifts to find true r
def A4map(x_in, y_in, z_in, t_in, L, Nt):
  r = np.sqrt((x_in**2 + y_in**2 + z_in**2 + t_in**2) * 0.8 \
              - 2.0 * (x_in * (y_in + z_in + t_in) \
                       + y_in * (z_in + t_in) + z_in * t_in) * 0.2)

  for x in [x_in + L, x_in, x_in - L]:
    xSq = x**2
    for y in [y_in + L, y_in, y_in - L]:
      xSq_ySq = xSq + y**2
      xpy = x + y
      xy = x * y
      for z in [z_in + L, z_in, z_in - L]:
        xSq_ySq_zSq = xSq_ySq + z**2
        xpypz = xpy + z
        xy_xz_yz = xy + x * z + y * z
        for t in [t_in + Nt, t_in, t_in - Nt]:
          test = np.sqrt((xSq_ySq_zSq + t**2) * 0.8 \
                         - 2.0 * (xy_xz_yz + xpypz * t) * 0.2)

          # Sanity check -- can be commented out for more speed
          x_a4 = (x - y) * invSq2
          y_a4 = (x + y - 2.0 * z) * invSq6
          z_a4 = (x + y + z - 3.0 * t) * invSq12
          t_a4 = (x + y + z + t) * invSq20
          check = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2 + t_a4**2)
          if np.fabs(test - check) > TOL:
            print "ERROR: %.6g isn't %.6g for (%d, %d, %d, %d)" \
                  % (check, test, x, y, z, t)
            sys.exit(1)

          # True r is the smallest
          # Try to avoid negative roundoff
          # so that we can see when periodic shifts are really necessary
          if test - r < -TOL:
            print "|(%d, %d, %d, %d)| = %.6g" % (x_in, y_in, z_in, t_in, r),
            print "--> |(%d, %d, %d, %d)| = %.6g" % (x, y, z, t, test)
            r = test
#  print "%d %d %d %d --> %.6g" % (x_in, y_in, z_in, t_in, r)
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
check = -1
for line in open(toOpen):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  # Format: hvy_pot: MAX_T = #, MAX_X = #
  # At least for now, forbid any component larger than MAX_X,
  # which is no larger than MAX_T for all cases we've run so far
  elif line.startswith('d_correlator_r: MAX_T = '):
    MAX_X = int((line.split())[6])

    # Find smallest r that is cut due to MAX_X
    MAX_r = 100.0 * MAX_X   # To be overwritten
    for y in range(MAX_X + 1):
      for z in range(MAX_X + 1):
        for t in range(MAX_X + 1):
          test = A4map(MAX_X + 1, y, z, t, L, Nt)
          if test < MAX_r:
            MAX_r = test
    print "MAX_r = %.6g" % MAX_r
    print "--------------------------------------"

  # Format: CORR_S x y z t dat
  elif line.startswith('CORR_S '):
    temp = line.split()
    x = int(temp[1])
    y = int(temp[2])
    z = int(temp[3])
    t = int(temp[4])
    dat = float(temp[5])
    if t > MAX_X:         # Drop large separations
      continue
    this_r = A4map(x, y, z, t, L, Nt)

    # Accumulate multiplicities if this_r < MAX_r
    # Try to avoid roundoff issues in comparison
    if this_r - MAX_r > -TOL:
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
  elif line.startswith('RUNNING COMPLETED'):
    check = 1
if check == -1:
  print toOpen, "did not complete"
  sys.exit(1)

# Sort by magnitude (column zero), not count
r = sorted(r, key=lambda x: x[0])
print "--------------------------------------"
for i in range(len(r)):   # Print r and average
  print "%.6g %.6g" % (r[i][0], r[i][2]  / r[i][1])
# ------------------------------------------------------------------

