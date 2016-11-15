#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Just combine Wilson loop data with the same |r|

# Parse argument: which files to consider (we'll do all of them)
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<tag>"
  sys.exit(1)
tag = str(sys.argv[1])
runtime = -time.time()

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

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
  r = 100.0 * L  # To be overwritten
  for x in [x_in + L, x_in, x_in - L]:
    for y in [y_in + L, y_in, y_in - L]:
      for z in [z_in + L, z_in, z_in - L]:
        test = np.sqrt((x**2 + y**2 + z**2) * 0.75 \
                       - (x * (y + z) + y * z) * 0.5)

        # Sanity check -- can be commented out for more speed
#        x_a4 = (x - y) * invSq2
#        y_a4 = (x + y - 2.0 * z) * invSq6
#        z_a4 = (x + y + z) * invSq12
#        check = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2)
#        if np.fabs(test - check) > TOL:
#          print "ERROR: %.4g isn't %.4g for (%d, %d, %d)" \
#                % (check, test, x, y, z)
#          sys.exit(1)

        # True r is the smallest
        if test < r:
          r = test
  return r
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through first output file to construct lookup table
# Also record list of values with corresponding multiplicity
# Decide where to cut r based on MAX_X
r = []          # List of two-component lists: first value, then count
files = 'Out/' + tag + '.*'
all_files = glob.glob(files)
for line in open(all_files[0]):
  if line.startswith('nx '):
    L = int((line.split())[1])
  # Format: hvy_pot: MAX_T = #, MAX_X = #
  elif line.startswith('hvy_pot: MAX_T '):
    temp = line.split()
    MAX_X = int(temp[6])
    num_y = 2 * MAX_X + 1
    lookup = np.empty((MAX_X + 1, num_y, num_y), dtype = np.float)

    # Decide where to cut r
    MAX_r = 100.0 * MAX_X   # To be overwritten
    for y in range(MAX_X + 1):
      for z in range(MAX_X + 1):
        test = A4map(MAX_X + 1, y, z, L)
        if test < MAX_r:
          MAX_r = test
#    print MAX_r

  # Format: *_LOOP n_x n_y n_z t dat
  elif '_LOOP ' in line:
    if ' 0 0 0 ' in line:     # Ignore degenerate (0, 0, 0)
      continue
    temp = line.split()
    if int(temp[4]) > 1:
      break                   # Only consider first t=1!
    x = int(temp[1])
    y = int(temp[2])
    z = int(temp[3])
    this_r = A4map(x, y, z, L)

    # Save this_r in lookup table so we can avoid calling A4map in the future
    lookup[x][y + MAX_X][z + MAX_X] = this_r

    # Accumulate multiplicities if this_r < MAX_r
    # Try to avoid roundoff issues in comparison
    if this_r - MAX_r > -TOL:
      continue
    done = -1
    for j in range(len(r)):
      if np.fabs(r[j][0] - this_r) < TOL:
        r[j][1] += 1
        done = 1
        break
    if done < 0:
      r.append([this_r, 1])
Npts = len(r)

# Sort by magnitude (column zero), not count
r = sorted(r, key=lambda x: x[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through all files matching input tag
ave = np.zeroes(Npts, dtype = np.float)
for filename in all_files:
  cfg = str(filename.split('.')[-1])       # Number after last .
  outfilename = 'Out/' + tag + '-naive.' + cfg
  outfile = open(outfilename, 'w')

  copyHeader = 1
  check = -1
  loop = "POT_LOOP"
  old_t = 1
  for line in open(filename):
    # Copy header to new output file
    if copyHeader > 0:
      print >> outfile, line.rstrip()

    if line.startswith('START '):   # Turn off copying after START line
      copyHeader = -1

    # Format: *_LOOP n_x n_y n_z t dat
    elif '_LOOP ' in line:
      if ' 0 0 0 ' in line:         # Ignore degenerate (0, 0, 0)
        continue
      temp = line.split()
      t = int(temp[4])

      # Average and reset if we've moved on to the next t
      # This will also take care of the three different types of loops
      if not t == old_t:
        for label in range(Npts):
          dat = ave[label] / r[label][1]
          print >> outfile, "%s %d %.4g %d %.6g" \
                            % (loop, label, r[label][0], old_t, dat)

          # Reset
          ave[label] = 0
        loop = str(temp[0])
        old_t = t

      # Accumulate data if this_r < MAX_r
      n_x = int(temp[1])
      n_y = int(temp[2])
      n_z = int(temp[3])
      this_r = lookup[x][y + MAX_X][z + MAX_X]
      if this_r - MAX_r > -TOL:
        continue
      done = -1
      for j in range(len(r)):
        if np.fabs(r[j][0] - this_r) < TOL:
          done = 1
          ave[j] += float(temp[5])
          break
      if done < 0:
        print "ERROR: displacement", this_r, "not found"
        sys.exit(1)

    elif line.startswith('RUNNING COMPLETED'):
      print >> outfile, line.rstrip()
      if check == 1:    # Check that we have one measurement per file
        print filename, "reports two measurements"
      check = 1

      # Average final data set and reset before moving on to the next file
      label = 0
      for a in range(num_x):
        for b in range(num_y):
          for c in range(num_z):
            if count[a][b][c] > 0:
              dat = ave[a][b][c] / count[a][b][c]
              print >> outfile, "%s %d %.4g %d %.6g" \
                                % (loop, label, rI[a][b][c], old_t, dat)
              label += 1

              # Reset
              ave[a][b][c] = 0.0
              count[a][b][c] = 0

  outfile.close()
  if check == -1:
    print filename, "did not complete"
    sys.exit(1)
# ------------------------------------------------------------------

