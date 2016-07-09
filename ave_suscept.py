#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA susceptibility averages
# with blocked standard errors

# For now only consider log-polar scalar field

# Currently the measurement output uses the naive |r| for the A4* lattice
# We just need to read, construct
#   \int dr r^{n + 3} \cO(x_0 + r) \cO(x_0) = \int dr r^{n + 3} C(r)
# for 0 <= n <= 3, and average

# Parse argument: the file to analyze
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
dist = []          # List of distances |r|
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
      dist.append(float(temp[3]))
  elif line.startswith('BLOCK_S 0 '):
    break             # Don't go through whole file yet

Npts = len(dist)

# Number of powers to consider: 0 <= n <= 4
Nsus = 4
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

    block = int(temp[1])
    label = int(temp[2])
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

# Integrated susceptibility estimates for all jk samples
susceptK = np.zeros((Nsus, Nblocks), dtype = np.float)
susceptS = np.zeros_like(susceptK)
for i in range(Nblocks):    # Jackknife samples
  # Need to sort data before we know which to integrate!
  K = np.zeros(Npts, dtype = columns)
  S = np.zeros_like(K)
  for x in range(Npts):
    K[x][0] = dist[x]
    K[x][1] = (Ktot[x] - Kdat[x][i]) / (Nblocks - 1.0)

    S[x][0] = dist[x]
    S[x][1] = (Stot[x] - Sdat[x][i]) / (Nblocks - 1.0)

  # Sort and record estimates -- note that r itself is not sorted
  # Include factor of dr since points are not evenly spaced!!!
  K = np.sort(K, order='r')
  S = np.sort(S, order='r')
  for n in range(Nsus):
    # !!! Skipping first three points to avoid discretization artifacts...
    for x in range(4, Npts):
      r = K[x][0]
      dr = K[x][0] - K[x - 1][0]
      susceptK[n][i] += np.power(r, n + 3) * K[x][1] * dr
      susceptS[n][i] += np.power(r, n + 3) * S[x][1] * dr
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing %d blocks of length %d MDTU" % (Nblocks, block_size)

outfile = open('results/suscept.dat', 'w')
print >> outfile, "# Analyzing %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# n Konishi err SUGRA err"

for n in range(Nsus):
  print >> outfile, "%d" % n,

  dat = np.array(susceptK[n])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.6g %.4g" % (ave, err),

  dat = np.array(susceptS[n])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.6g %.4g" % (ave, err)
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
