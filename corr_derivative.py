#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Compute the numerical derivative of the radial Konishi and SUGRA correlators
# with blocked standard errors

# For now only consider log-polar operator
# Subtraction should vanish in the derivative -- possible consistency check

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us which files to analyze
#   (for example, "corr", "konishi0", "subtracted0")
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct list of which configurations have been analyzed
cfgs = []
files = 'Out/' + tag + '.*'
for filename in glob.glob(files):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs and cfg > cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files", files, "found"
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]

# Cycle through first output file to associate scalar distances
# with the corresponding label in the output
r = []          # List of r
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  # Format: CORR_K label r tag1 tag2 dat_vev dat_vol
  # Format: CORR_K label r dat_vev dat_vol
  if line.startswith('CORR_K '):
    temp = line.split()
    tag1 = int(temp[3])
    tag2 = int(temp[4])
    if tag1 == 0 and tag2 == 0:
      tr = float(temp[2])
      if tr == 0:     # Skip r=0!
        continue
      r.append(tr)
  elif line.startswith('CORR_S '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
# For now only consider log-polar operator with ensemble subtraction
Kdat = [[] for x in range(Npts)]
Sdat = [[] for x in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [[] for x in range(Npts)]
tS = [[] for x in range(Npts)]
for i in range(Npts):
  tK[i] = 0.0
  tS[i] = 0.0

for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for i in range(Npts):
      Kdat[i].append(tK[i] / float(count))
      Sdat[i].append(tS[i] / float(count))
      tK[i] = 0.0
      tS[i] = 0.0
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  for line in open(toOpen[0]):
    # Format: CORR_? label r tag1 tag2 dat_vev dat_vol
    if line.startswith('CORR'):
      if line.startswith('CORR_K 0 0 0 0 '):
        count += 1            # Only increment once per measurement!

      temp = line.split()
      if float(temp[2]) == 0: # Skip r=0!
        continue

      # Shift label on r since we skip r=0...
      label = int(temp[1]) - 1
      tag1 = int(temp[3])
      tag2 = int(temp[4])
      dat = float(temp[5])

      # For now only consider log-polar operator with ensemble subtraction
      if line.startswith('CORR_K ') and tag1 == 0 and tag2 == 0:
        tK[label] += dat
      elif line.startswith('CORR_S ') and tag1 == 0 and tag2 == 0:
        tS[label] += dat

    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print infile, "reports two measurements"
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for i in range(Npts):
    Kdat[i].append(tK[i] / float(count))
    Sdat[i].append(tS[i] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
# Require multiple blocks for now
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

columns = [('r', float), ('dat', float)]
Ktot   = np.array([sum(Kdat[x]) for x in range(Npts)])
Stot   = np.array([sum(Sdat[x]) for x in range(Npts)])

# Effective dimension estimates for all jk samples
derivK = np.empty((Npts - 1, Nblocks), dtype = np.float)
derivS = np.empty((Npts - 1, Nblocks), dtype = np.float)
for i in range(Nblocks):    # Jackknife samples
  # Need to sort data before we know which to divide!
  K = np.zeros(Npts, dtype = columns)
  S = np.zeros_like(K)
  for x in range(Npts):
    K[x][0] = r[x]
    K[x][1] = (Ktot[x] - Kdat[x][i]) / (Nblocks - 1.0)

    S[x][0] = r[x]
    S[x][1] = (Stot[x] - Sdat[x][i]) / (Nblocks - 1.0)

  # Sort and record estimates -- note that r itself is not sorted
  K = np.sort(K, order='r')
  S = np.sort(S, order='r')
  for x in range(Npts - 1):
    derivK[x][i] = (K[x + 1][1] - K[x][1]) / (K[x + 1][0] - K[x][0])
    derivS[x][i] = (S[x + 1][1] - S[x][1]) / (S[x + 1][0] - S[x][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing %d blocks of length %d MDTU" % (Nblocks, block_size)
Kfile = open('results/deriv_konishi.dat', 'w')
print >> Kfile, "# Analyzing %d blocks of length %d MDTU" \
                % (Nblocks, block_size)

Sfile = open('results/deriv_sugra.dat', 'w')
print >> Sfile, "# Analyzing %d blocks of length %d MDTU" \
                % (Nblocks, block_size)

for x in range(Npts - 1):
  mid_r = 0.5 * (K[x][0] + K[x + 1][0])

  # Konishi
  ave = np.average(derivK[x])
  var = (Nblocks - 1.0) * np.sum((derivK[x] - ave)**2) / float(Nblocks)
  print >> Kfile, "%.4g %.6g %.4g" % (mid_r, -0.5 * ave, 0.5 * np.sqrt(var))

  # SUGRA
  ave = np.average(derivS[x])
  var = (Nblocks - 1.0) * np.sum((derivS[x] - ave)**2) / float(Nblocks)
  print >> Sfile, "%.4g %.6g %.4g" % (mid_r, -0.5 * ave, 0.5 * np.sqrt(var))
Kfile.close()
Sfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
