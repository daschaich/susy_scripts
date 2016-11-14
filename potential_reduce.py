#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Just combine Wilson loop data with the same r_I (unweighted averages)
# r_I computed with vegas as described in notes, and copied into a lookup table

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

# Lookup table is customized for lattice volume 8nt24,
# but other volumes should work up to potentially missing large-r points
temp = os.getcwd()
L = -1
if '_8nt24' in temp:
  L = 8
elif '_12nt24' in temp:
  L = 12
elif '_16nt24' in temp:
  L = 16
else:
  print "Error: Unrecognized lattice volume in", temp
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# The lookup table, plus tables of the same size for running averages
if L == 8:
  num_x = 3;    num_y = 3;      num_z = 4
elif L == 12:
elif L == 16:
count = np.zeros((num_x, num_y, num_z), dtype = np.int)
ave = np.zeros((num_x, num_y, num_z), dtype = np.float)
rI = np.empty((num_x, num_y, num_z), dtype = np.float)
for a in range(num_x):
  for b in range(num_y):
    for c in range(num_z):      # Negative default will make it obvious
      rI[a][b][c] = -1.0        # if we overlooked a displacement
rI[0][0][1] = 0.919
rI[0][1][1] = 1.411
rI[0][0][2] = 1.752
rI[1][1][1] = 1.810
rI[0][1][2] = 2.154
rI[1][1][2] = 2.461
rI[0][0][3] = 2.776
rI[0][2][2] = 2.792
rI[1][2][2] = 3.022
rI[0][1][3] = 3.051
rI[1][1][3] = 3.272
rI[2][2][2] = 3.504
rI[0][2][3] = 3.568
rI[1][2][3] = 3.739
if L > 10:
  rI[0][0][4] = 3.836
  rI[0][1][4] = 4.017
  rI[2][2][3] = 4.147
  rI[1][1][4] = 4.175
  rI[0][3][3] = 4.227
  rI[1][3][3] = 4.364
  rI[0][2][4] = 4.425
  rI[1][2][4] = 4.557
  rI[2][3][3] = 4.716
  rI[0][0][5] = 4.879
  rI[2][2][4] = 4.900
  rI[0][3][4] = 4.984
  rI[0][1][5] = 5.006
  rI[1][3][4] = 5.096
  rI[1][1][5] = 5.125
  rI[3][3][3] = 5.226
  rI[0][2][5] = 5.333
  rI[2][3][4] = 5.399
  rI[1][2][5] = 5.439
  rI[0][4][4] = 5.651
  rI[2][2][5] = 5.728
  rI[1][4][4] = 5.746
  rI[0][3][5] = 5.809
  rI[3][3][4] = 5.854
  rI[1][3][5] = 5.902
if L > 15:
  rI[0][0][6] = 5.906
  rI[0][1][6] = 6.003
  rI[2][4][4] = 6.014
  rI[1][1][6] = 6.097
  rI[2][3][5] = 6.165
  rI[0][2][6] = 6.272
  rI[1][2][6] = 6.359
  rI[0][4][5] = 6.396
  rI[3][4][4] = 6.425
  rI[1][4][5] = 6.478
  rI[3][3][5] = 6.569
  rI[2][2][6] = 6.606
  rI[0][3][6] = 6.680
  rI[2][4][5] = 6.716
  rI[1][3][6] = 6.760
  rI[0][0][7] = 6.922
  rI[4][4][4] = 6.952
  rI[2][3][6] = 6.989
  rI[0][1][7] = 7.002
  rI[0][5][5] = 7.068
  rI[1][1][7] = 7.079
  rI[3][4][5] = 7.087
  rI[1][5][5] = 7.142
  rI[0][4][6] = 7.199
  rI[0][2][7] = 7.229
  rI[1][4][6] = 7.272
  rI[1][2][7] = 7.303
  rI[3][3][6] = 7.349
  rI[2][5][5] = 7.356
  rI[2][4][6] = 7.483
  rI[2][2][7] = 7.517
  rI[4][4][5] = 7.569
  rI[0][3][7] = 7.584
  rI[1][3][7] = 7.6563
  rI[3][5][5] = 7.696
  rI[0][5][6] = 7.806
  rI[3][4][6] = 7.818
  rI[2][3][7] = 7.856
  rI[1][5][6] = 7.872
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

