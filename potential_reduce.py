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

# Only constructed lookup table for volumes up to 8nt24 so far!!!
temp = os.getcwd()
if not '_8nt24' in temp:
  print "ERROR: Can only handle volumes up to 8nt24"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# The lookup table, plus tables of the same size for running averages
# We can ignore the sign of r1, but not r2 or r3
# r2 ranges from -6 through 8: total 15, shifted by 6
# r3 ranges from -3 through 9: total 13, shifted by 3
# awk '{print "rI["$1"]["$2+6"]["$3+3"] = "$5}' <file>
max_x = 5;    max_y = 15;     max_z = 13
count = np.zeros((max_x, max_y, max_z), dtype = np.int)
ave = np.zeros((max_x, max_y, max_z), dtype = np.float)
rI = np.empty((max_x, max_y, max_z), dtype = np.float)
for a in range(max_x):
  for b in range(max_y):
    for c in range(max_z):      # Negative default will make it obvious
      rI[a][b][c] = -1.0        # if we overlooked a vector below
rI[0][6][6] = 0.542
rI[0][4][4] = 0.790
rI[1][7][4] = 0.802
rI[0][8][5] = 0.823
rI[1][5][5] = 0.824
rI[0][6][9] = 1.052
rI[2][6][3] = 1.276
rI[0][4][7] = 1.303
rI[1][5][8] = 1.361
rI[0][8][8] = 1.368
rI[1][7][7] = 1.368
rI[1][9][3] = 1.417
rI[0][10][4] = 1.484
rI[2][6][6] = 1.625
rI[2][8][2] = 1.713
rI[2][4][4] = 1.713
rI[0][2][5] = 1.795
rI[1][9][6] = 1.801
rI[0][10][7] = 1.811
rI[1][3][6] = 1.853
rI[3][7][1] = 1.976
rI[3][5][5] = 1.977
rI[3][7][4] = 1.991
rI[3][5][2] = 1.993
rI[1][11][2] = 2.027
rI[2][10][4] = 2.084
rI[2][8][5] = 2.085
rI[1][11][5] = 2.091
rI[2][10][1] = 2.105
rI[2][2][5] = 2.105
rI[2][4][7] = 2.108
rI[0][6][12] = 2.223
rI[0][12][3] = 2.345
rI[3][9][3] = 2.388
rI[3][3][3] = 2.388
rI[0][8][11] = 2.415
rI[1][5][11] = 2.425
rI[0][10][10] = 2.468
rI[1][1][7] = 2.473
rI[3][5][8] = 2.497
rI[0][12][6] = 2.524
rI[3][7][7] = 2.541
rI[3][9][6] = 2.562
rI[2][10][7] = 2.584
rI[0][0][6] = 2.615
rI[3][3][6] = 2.633
rI[3][9][0] = 2.634
rI[2][2][8] = 2.669
rI[4][6][3] = 2.678
rI[2][6][9] = 2.704
rI[0][2][8] = 2.733
rI[1][11][8] = 2.751
rI[2][4][10] = 2.815
rI[2][12][3] = 2.830
rI[1][13][4] = 2.866
rI[3][11][2] = 2.973
rI[3][1][4] = 2.975
rI[4][8][2] = 2.981
rI[4][4][4] = 2.982
rI[1][9][9] = 2.982
rI[4][6][6] = 3.003
rI[4][6][0] = 3.004
rI[0][4][10] = 3.009
rI[0][12][9] = 3.079
rI[1][13][1] = 3.104
rI[2][12][6] = 3.127
rI[0][14][5] = 3.135
rI[2][8][8] = 3.133
rI[2][0][6] = 3.133
rI[2][12][0] = 3.135
rI[1][13][7] = 3.138
rI[0][14][2] = 3.138
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
        for a in range(max_x):
          for b in range(max_y):
            for c in range(max_z):
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

      # Convert (n_x, n_y, n_z) to (r_1, r_2, r_3) tag,
      # shifting r_2 and r_3 to keep them non-negative
      # Instead of computing too-large distances,
      # use rI array to figure out which are cut
      n_x = int(temp[1])
      n_y = int(temp[2])
      n_z = int(temp[3])
      rtag = [np.abs(n_x - n_y),
              6 + n_x + n_y - 2 * n_z,
              3 + n_x + n_y + n_z]
      if rtag[0] >= max_x or rtag[1] >= max_y or rtag[2] >= max_z \
                          or rI[rtag[0]][rtag[1]][rtag[2]] < 0:
        continue

      ave[rtag[0]][rtag[1]][rtag[2]] += float(temp[5])
      count[rtag[0]][rtag[1]][rtag[2]] += 1

    elif line.startswith('RUNNING COMPLETED'):
      print >> outfile, line.rstrip()
      if check == 1:    # Check that we have one measurement per file
        print filename, "reports two measurements"
      check = 1

      # Average final data set and reset before moving on to the next file
      label = 0
      for a in range(max_x):
        for b in range(max_y):
          for c in range(max_z):
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

