#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compute SUGRA correlator effective masses, using a blocked jackknife
# For C(t) = A (e^{-Mt} + e^{-M(T-t)}),
# arccosh M = [C(t-1) + C(t+1)] / 2C(t)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us whether to analyze "corr" or "stout" files
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
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

# Extract Nt from first output file
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nt'):
    Nt = int((line.split())[1])
    break

# We exploit t <--> Nt - t symmetry
tmax = int(Nt / 2)   # Assume Nt is even
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# For now just average over all 25 SUGRA components
# This constant factor will cancel out in the effective mass ratio
# We lose two points in the ratio [C(t-1) + C(t+1)] / 2C(t)
effm = [[] for x in range(tmax - 1)]
warnings = 0

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tM = [0 for x in range(tmax + 1)]   # Include x=tmax
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for t in range(1, tmax):
      temp = (tM[t - 1] + tM[t + 1]) / (2.0 * tM[t])
      if temp < 1:    # Skip any malformed effective masses
#        print "WARNING: Malformed m_eff for t=%d at block %d" \
#              % (t, len(effm[t - 1]) + 1)
        warnings += 1
      else:
        effm[t - 1].append(np.arccosh(temp))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tM = [0 for x in range(tmax + 1)]   # Include x=tmax

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: SUGRA t dat
    if line.startswith('SUGRA '):
      temp = line.split()
      time = int(temp[1])
      dat = float(temp[2])
      tM[time] += dat

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for t in range(1, tmax):
    temp = (tM[t - 1] + tM[t + 1]) / (2.0 * tM[t])
    if temp < 1:    # Skip any malformed effective masses
#      print "WARNING: Malformed m_eff for t=%d at block %d" \
#            % (t, len(effm[t - 1]) + 1)
      warnings += 1
    else:
      effm[t - 1].append(np.arccosh(temp))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(block_data[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error
# Allow different numbers of measurements for each t
print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
print "%d of %d eff_m results were malformed" % (warnings, Nblocks * len(effm))
outfile = open('results/sugra_effm.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# %d of %d eff_m results were malformed" \
                  % (warnings, Nblocks * len(effm))
for t in range(1, tmax):
  Nblocks = len(effm[t - 1])
  if Nblocks <= 1:
    print "ERROR: not enough valid measurements for t =", t
    continue

  dat = np.array(effm[t - 1])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> outfile, "%d %.6g %.4g # %d" % (t, ave, err, Nblocks)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
