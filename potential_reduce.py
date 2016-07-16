#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Just combine Wilson loop data with the same r_I (unweighted averages)
# r_I computed with vegas as described in notes,
# and copied into a lookup table indexed by the coefficients of
# 1/sqrt(2), 1/sqrt(6) and 1/sqrt(12) in the naive r on an A4* timeslice

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

# Lookup table is customized for lattice volume 8nt24
temp = os.getcwd()
if not '_8nt24' in temp:
  print "ERROR: Customized for lattice volume 8nt24"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# The lookup table, plus tables of the same size for running averages
num_x = 3;    num_y = 3;      num_z = 4
count = np.zeros((num_x, num_y, num_z), dtype = np.int)
ave = np.zeros((num_x, num_y, num_z), dtype = np.float)
rI = np.empty((num_x, num_y, num_z), dtype = np.float)
for a in range(num_x):
  for b in range(num_y):
    for c in range(num_z):      # Negative default will make it obvious
      rI[a][b][c] = -1.0        # if we overlooked a vector below
rI[0][0][1] = 0.925
rI[0][1][1] = 1.443
rI[1][1][1] = 1.827
rI[0][0][2] = 1.856
rI[0][1][2] = 2.215
rI[1][1][2] = 2.490
rI[0][2][2] = 2.838
rI[0][0][3] = 2.890
rI[1][2][2] = 3.042
rI[0][1][3] = 3.117
rI[1][1][3] = 3.309
rI[2][2][2] = 3.514
rI[0][2][3] = 3.605
rI[1][2][3] = 3.761
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through all files matching input tag
files = 'Out/' + tag + '.*'
for filename in glob.glob(files):
  cfg = str(filename.split('.')[-1])       # Number after last .
  outfilename = 'Out/' + tag + '-reduced.' + cfg
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
      if ' 0 0 0 ' in line:    # Ignore degenerate (0, 0, 0)
        continue
      temp = line.split()
      t = int(temp[4])

      # Average and reset if we've moved on to the next t
      # This will also take care of the three different types of loops
      if not t == old_t:
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
        loop = str(temp[0])
        old_t = t

      # Convert (n_x, n_y, n_z) displacements to canonical form
      # (absolute values in increasing order)
      n_x = np.abs(int(temp[1]))
      n_y = np.abs(int(temp[2]))
      n_z = np.abs(int(temp[3]))
      n = np.sort([n_x, n_y, n_z])
      if n[0] >= num_x or n[1] >= num_y or n[2] >= num_z \
                       or rI[n[0]][n[1]][n[2]] < 0:
        continue

      ave[n[0]][n[1]][n[2]] += float(temp[5])
      count[n[0]][n[1]][n[2]] += 1

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

