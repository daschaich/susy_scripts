#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Extract the blocked Konishi and SUGRA effective dimension
#   De_Eff = (1/2) log[C(r1) / C(r2)] / log[r2 / r1]
# with blocked standard errors

# For now only consider log-polar scalar field

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

# Take a first pass through the file to read the number of blocks
# and associate scalar distances with the corresponding label in the output
r = []          # List of r
for line in open(toOpen):
  if line.startswith('Nblock '):
    Nblocks = int((line.split())[1])

  # !!! Assume measurements always separated by 10 MDTU
  elif line.startswith('Nmeas '):
    block_size = 10 * int((line.split())[1])

  # Format: BLOCK_K block# label r tag1 tag2 dat
  elif line.startswith('BLOCK_K 0 '):
    temp = line.split()
    tag1 = int(temp[4])
    tag2 = int(temp[5])
    if tag1 == 0 and tag2 == 0:
      tr = float(temp[3])
      if tr == 0:     # Skip r=0!
        continue
      r.append(tr)
  elif line.startswith('BLOCK_S 0 '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read in blocked measurements of each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
# For now only consider log-polar scalar field
Kdat = np.empty((Npts, Nblocks), dtype = np.float)
Sdat = np.empty_like(Kdat)

check = -1
for line in open(toOpen):
  # Format: BLOCK_K block# label r tag1 tag2 dat
  if line.startswith('BLOCK'):
    temp = line.split()
    if float(temp[3]) == 0: # Skip r=0!
      continue

    # Shift label on r since we skip r=0...
    block = int(temp[1])
    label = int(temp[2]) - 1
    tag1 = int(temp[4])
    tag2 = int(temp[5])
    if not tag1 == 0 or not tag2 == 0:
      continue
    dat = float(temp[6])

    if line.startswith('BLOCK_K '):
      Kdat[label][block] = dat
    elif line.startswith('BLOCK_S '):
      Sdat[label][block] = dat

  elif line.startswith('RUNNING COMPLETED'):
    if check == 1:    # Check that we have one measurement per file
      print infile, "reports two measurements"
    check = 1
if check == -1:
  print toOpen, "did not complete"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
# Require multiple blocks for now
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

columns = [('r', float), ('dat', float)]
Ktot    = np.array([sum(Kdat[x]) for x in range(Npts)])
Stot    = np.array([sum(Sdat[x]) for x in range(Npts)])

# Effective dimension estimates for all jk samples
DeK = np.empty((Npts - 1, Nblocks), dtype = np.float)
DeS = np.empty_like(DeK)
for i in range(Nblocks):    # Jackknife samples
  # Need to sort data before we know which to divide!
  K = np.zeros(Npts, dtype = columns)
  S = np.zeros_like(K)
  for x in range(Npts):
    K[x][0] = r[x]
    K[x][1] = (Ktot[x] - Kdat[x][i]) / (Nblocks - 1.0)

    S[x][0] = r[x]
    S[x][1] = (Stot[x] - Sdat[x][i]) / (Nblocks - 1.0)

    if K[x][1] < 0:
      print "Warning: Negative C_K(%.4g) = %.4g" % (K[x][0], K[x][1])
    if S[x][1] < 0:
      print "Warning: Negative C_S(%.4g) = %.4g" % (S[x][0], S[x][1])

  # Sort and record estimates -- note that r itself is not sorted
  K = np.sort(K, order='r')
  S = np.sort(S, order='r')
  for x in range(Npts - 1):
    DeK[x][i] = np.log(K[x][1] / K[x + 1][1])
    DeK[x][i] *= 0.5 / np.log(K[x + 1][0] / K[x][0])
    DeS[x][i] = np.log(S[x][1] / S[x + 1][1])
    DeS[x][i] *= 0.5 / np.log(S[x + 1][0] / S[x][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing %d blocks of length %d MDTU" % (Nblocks, block_size)

outfile = open('results/eff_delta.dat', 'w')
print >> outfile, "# Analyzing %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# r    De_K      err     De_S      err"
for x in range(Npts - 1):
  mid_r = 0.5 * (K[x][0] + K[x + 1][0])

  # Konishi
  ave = np.average(DeK[x])
  var = (Nblocks - 1.0) * np.sum((DeK[x] - ave)**2) / float(Nblocks)
  print >> outfile, "%.4g %.6g %.4g" % (mid_r, ave, np.sqrt(var)),

  # SUGRA
  ave = np.average(DeS[x])
  var = (Nblocks - 1.0) * np.sum((DeS[x] - ave)**2) / float(Nblocks)
  print >> outfile, "%.6g %.4g" % (ave, np.sqrt(var))
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
