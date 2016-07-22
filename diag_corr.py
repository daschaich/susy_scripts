#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Run jackknifed diagonalization of radial correlators

# For now only consider log-polar scalar field
# For now hard-code 2x2 matrices
Nops = 2

# Parse argument: which file to analyze
# This file already handles the thermalization cut and blocking
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Take a first pass through the file to read the number of blocks,
# count how many scalar distances we will consider, and record their order
r = []          # List of r
for line in open(toOpen):
  if line.startswith('Nblock '):
    Nblocks = int((line.split())[1])

  # !!! Assume measurements always separated by 10 MDTU
  elif line.startswith('Nmeas '):
    block_size = 10 * int((line.split())[1])

  # Format: BLOCK block# label r tag1 tag2 dat
  elif line.startswith('BLOCK 0 '):
    temp = line.split()
    tag1 = int(temp[4])
    tag2 = int(temp[5])
    if tag1 == 0 and tag2 == 0:   # Only count each r once
      r.append(float(temp[3]))
  elif line.startswith('CORR 0 '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read in blocked measurements of each correlator
# For now only consider log-polar scalar field
dat = np.empty((Npts, Nops, Nblocks), dtype = np.float)

check = -1
mat = np.zeros((Nops, Nops), dtype = np.float)
oldBlock = 0
rCount = 0
for line in open(toOpen):
  # Format: BLOCK block# label r tag1 tag2 dat
  if line.startswith('BLOCK '):
    temp = line.split()
    block = int(temp[1])
    tag1 = int(temp[4])
    tag2 = int(temp[5])

    # Check and reset counters if we've moved on to the next block
    if not block == oldBlock:
      if not rCount == Npts:
        print "Wrong number of points in block", oldBlock
      oldBlock = block
      rCount = 0

    mat[tag1][tag2] = float(temp[6])

    # If we're done with this r for this block, diagonalize
    # Check for complex eigenvalues (negative eigenvalues should be okay)
    if tag1 == Nops - 1 and tag2 == Nops - 1:
      # Symmetrize by hand
      tr = 0.5 * (mat[0][1] + mat[1][0])
      mat[0][1] = tr
      mat[1][0] = tr
      eig, vecs = np.linalg.eig(mat)

      re_eig = np.zeros(Nops, dtype = np.float)
      for i in range(Nops):
        if np.absolute(eig[i].imag) > 1e-8:
          print "ERROR: Complex eigenvalue"
          print eig
#          sys.exit(1)
        re_eig[i] = eig[i].real
      temp = np.sort(re_eig)       # Orders from smallest to largest
      re_eig = temp[::-1]          # Reverse order
      dat[rCount][0][block] = re_eig[0]
      dat[rCount][1][block] = re_eig[1]
      mat = np.zeros((Nops, Nops), dtype = np.float)    # Reset

      # Absolute values still lead to error bars consistent with zero
#      dat[rCount][0][block] = np.amax(np.absolute(eig))
#      dat[rCount][1][block] = np.amin(np.absolute(eig))
      rCount += 1

  elif line.startswith('RUNNING COMPLETED'):
    if check == 1:    # Check that we have one measurement per file
      print infile, "reports two measurements"
    check = 1
if check == -1:
  print toOpen, "did not complete"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)
print "Diagonalizing correlators built from %d operators" % Nops

outfile = open('results/corr.dat', 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# Diagonalizing correlators built from %d operators" % Nops
print >> outfile, "# r    Konishi   err       SUGRA     err"

for x in range(Npts):
  ave = np.mean(dat[x][0], dtype = np.float64)
  err = np.std(dat[x][0], dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.4g %.6g %.4g" % (r[x], ave, err),

  ave = np.mean(dat[x][1], dtype = np.float64)
  err = np.std(dat[x][1], dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.6g %.4g" % (ave, err)
outfile.close()

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime
# ------------------------------------------------------------------
