#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Jackknife solutions of generalized eigenvalue problem (GEVP)
# Given radial correlators C(r), find eigenvalues of
#   C(x0)^{-1} C(r) ~ (r / x0)^{-Delta}

# For now only consider log-polar scalar field
# For now hard-code 2x2 matrices
Nops = 2
# For now just print all lambda(r, x0) rather than fitting them

# Parse argument: which file to analyze
# This file already handles the thermalization cut and blocking
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Diagonalize all C(x0)^{-1} C(x) for x > x0
def diagonalize(mats, la_dat, Npts, block):
  for x0 in range(Npts - 1):
    Cinv = np.linalg.inv(mats[x0])
    for x in range(x0 + 1, Npts):
      toDiag = np.dot(Cinv, mats[x])
      print Cinv
      print mats[x]
      print toDiag
      eig, vecs = np.linalg.eig(toDiag)

      re_eig = np.empty(Nops, dtype = np.float)
      for i in range(Nops):
        if np.absolute(eig[i].imag) > 1e-8:
          print "ERROR: Complex eigenvalue"
          print eig
#          sys.exit(1)
        # Absolute value leads to non-zero asymptotic value
        re_eig[i] = eig[i].real

      temp = np.sort(re_eig)       # Orders from smallest to largest
      re_eig = temp[::-1]          # Reverse order
      for i in range(Nops):
        la_dat[x0][x][i][block] = re_eig[i]
  sys.exit(0)
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
  elif line.startswith('BLOCK 1 '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read in blocked measurements of each correlator
# For now only consider log-polar scalar field
la_dat = np.empty((Npts, Npts, Nops, Nblocks), dtype = np.float)

check = -1
mats = np.empty((Npts, Nops, Nops), dtype = np.float)
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

    mats[rCount][tag1][tag2] = float(temp[6])
    # !!!For now symmetrize by hand...
    if tag1 == Nops - 1 and tag2 == Nops - 1:
#      tr = 0.5 * (mats[rCount][0][1] + mats[rCount][1][0])
#      mats[rCount][0][1] = tr
#      mats[rCount][1][0] = tr
      rCount += 1

    # If we're done with all r for this block
    # Compute and diagonalize C(x0)^{-1} C(r) for r > x0
    # Check for complex eigenvalues (negative eigenvalues should be okay)
    if rCount == Npts - 1 and tag1 == Nops - 1 and tag2 == Nops - 1:
      diagonalize(mats, la_dat, Npts, block)
      mats = np.empty((Npts, Nops, Nops), dtype = np.float)    # Reset

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

outfile = open('results/GEVP.dat', 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# Diagonalizing correlators built from %d operators" % Nops
print >> outfile, "# r    r0    Konishi   err       SUGRA     err"

for x0 in range(Npts - 1):
  for x in range(x0 + 1, Npts):
    ave = np.mean(dat[x][x0][0], dtype = np.float64)
    err = np.std(dat[x][x0][0], dtype = np.float64) / np.sqrt(Nblocks - 1.0)
    print >> outfile, "%.4g %.4g %.6g %.4g" % (r[x0], r[x], ave, err),

    ave = np.mean(dat[x][x0][1], dtype = np.float64)
    err = np.std(dat[x][x0][1], dtype = np.float64) / np.sqrt(Nblocks - 1.0)
    print >> outfile, "%.6g %.4g" % (ave, err)
outfile.close()

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime
# ------------------------------------------------------------------
