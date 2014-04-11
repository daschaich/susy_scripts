#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked plaquette and R-symmetry-modified Wilson loop averages
# with blocked standard errors

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct list of which configurations have been analyzed
cfgs = []
for filename in glob.glob('Out/corr.*'):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs and cfg > cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files Out/corr.* found"
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each observable
# P = usual plaquette, F = modified loop
# D = difference (F - P), R = ratio (F - P)/(F + P)
Pdat = []
Fdat = []
Ddat = []
Rdat = []

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tP = 0.0
tF = 0.0
tD = 0.0
tR = 0.0
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    Pdat.append(tP / float(count))
    Fdat.append(tF / float(count))
    Ddat.append(tD / float(count))
    Rdat.append(tR / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tP = 0.0
    tF = 0.0
    tD = 0.0
    tR = 0.0

  # Running averages
  dat = 0.0
  filename = 'Out/corr.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: GMES Re(Poly) Im(Poly) cg_iters ss_plaq st_plaq S_B
    # Average over space--space and space--time plaquettes
    # Until 17 January 2014, the plaquettes didn't involve the diagonal link
    if line.startswith('GMES '):
      temp = line.split()
      sav = (float(temp[4]) + float(temp[5])) / 2.0
      tP += sav
      count += 1    # Only increment once per measurement!
    # Format: FLAVOR a b dat
    elif line.startswith('FLAVOR '):
      temp = line.split()
      dat += float(temp[3])
      # Average before comparing to plaquette, after this last one
      if int(temp[1]) == 4 and int(temp[2]) == 3:
        dat += float(temp[3])
        dat /= 20.0
        tF += dat
        tD += (dat - sav)
        tR += (dat - sav) / (dat + sav)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  Pdat.append(tP / float(count))
  Fdat.append(tF / float(count))
  Ddat.append(tD / float(count))
  Rdat.append(tR / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Pdat)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
outfile = open('results/flavor.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
print >> outfile, "# diff err rel err flavor err plaq err"

dat = np.array(Ddat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Rdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Fdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Pdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g # %d" % (ave, err, Nblocks)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
