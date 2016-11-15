#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Just combine Wilson loop data with the same r_I (unweighted averages)
# r_I computed with 8nt24 discrete Fourier transform,
# and copied into a lookup table

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

# Lookup table is specific to lattice volume
# Other volumes need to be re-run
temp = os.getcwd()
if '_8nt24' in temp:
  vol = '8nt24'
  num_x = 2;    num_y = 3;      num_z = 4
elif '_12nt24' in temp:
  vol = '12nt24'
  num_x = 3;    num_y = 4;      num_z = 6
else:
  print "ERROR: Lattice volume not yet implemented"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# The lookup table, plus tables of the same size for running averages
count = np.zeros((num_x, num_y, num_z), dtype = np.int)
ave = np.zeros((num_x, num_y, num_z), dtype = np.float)
rI = np.empty((num_x, num_y, num_z), dtype = np.float)
for a in range(num_x):
  for b in range(num_y):
    for c in range(num_z):      # Negative default will make it obvious
      rI[a][b][c] = -1.0        # if we overlooked a displacement

if vol == '8nt24':
  rI[0][0][1] = 0.913
  rI[0][1][1] = 1.385
  rI[0][0][2] = 1.694
  rI[1][1][1] = 1.752
  rI[0][1][2] = 2.045
  rI[1][1][2] = 2.298
  rI[0][0][3] = 2.468
  rI[0][2][2] = 2.541
  rI[0][1][3] = 2.660
  rI[1][2][2] = 2.707

elif vol == '12nt24':
  rI[0][0][1] =  0.926
  rI[0][1][1] =  1.434
  rI[0][0][2] =  1.794
  rI[1][1][1] =  1.856
  rI[0][1][2] =  2.231
  rI[1][1][2] =  2.576
  rI[0][0][3] =  2.928
  rI[0][2][2] =  2.956
  rI[1][2][2] =  3.230
  rI[0][1][3] =  3.254
  rI[1][1][3] =  3.523
  rI[2][2][2] =  3.824
  rI[0][2][3] =  3.888
  rI[1][2][3] =  4.110
  rI[0][0][4] =  4.165
  rI[0][1][4] =  4.398
  rI[1][1][4] =  4.606
  rI[2][2][3] =  4.649
  rI[0][3][3] =  4.736
  rI[1][3][3] =  4.926
  rI[0][2][4] =  4.939
  rI[1][2][4] =  5.123
  rI[0][0][5] =  5.247
  rI[0][1][5] =  5.420
  rI[2][3][3] =  5.424
  rI[1][1][5] =  5.585
  rI[2][2][4] =  5.609
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through all files matching input tag
files = 'Out/' + tag + '.*'
for filename in glob.glob(files):
  cfg = str(filename.split('.')[-1])       # Number after last .
  outfilename = 'Out/' + tag + '-FVreduced.' + cfg
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

  print >> outfile, "RUNNING COMPLETED"
  outfile.close()
  if check == -1:
    print filename, "did not complete"
    sys.exit(1)
# ------------------------------------------------------------------

