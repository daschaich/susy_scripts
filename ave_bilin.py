#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked fermion bilinear and susy transformation averages
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
  if cfg not in cfgs and cfg >= cut:
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
# Ignore stochastic noise in measurements
# S = SU(2) part of bilinear, U = U(1) part of bilinear, G = gauge piece,
# O = susy transform (G - S - U), R = relative (G - S - U)/(G + S + U)
Sdat = []
Udat = []
Gdat = []
Odat = []
Rdat = []

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tS = 0.0
tU = 0.0
tG = 0.0
tO = 0.0
tR = 0.0
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    Sdat.append(tS / float(count))
    Udat.append(tU / float(count))
    Gdat.append(tG / float(count))
    Odat.append(tO / float(count))
    Rdat.append(tR / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tS = 0.0
    tU = 0.0
    tG = 0.0
    tO = 0.0
    tR = 0.0

  # Running averages
  filename = 'Out/corr.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: BILIN dat imag (last should average to zero)
    if line.startswith('BILIN '):
      temp = line.split()
      count += 1    # Only increment once per measurement!
      dat = float(temp[1])
      tS += dat
    # Format: SUSY f_dat imag g_dat diff
    # As above, imaginary component should average to zero
    elif line.startswith('SUSY '):
      temp = line.split()
      tU += float(temp[1]) - dat    # Only accumulate U(1) piece
      tG += float(temp[3])
      tO += float(temp[4])
      tR += float(temp[4]) / (float(temp[1]) + float(temp[3]))

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  Sdat.append(tS / float(count))
  Udat.append(tU / float(count))
  Gdat.append(tG / float(count))
  Odat.append(tO / float(count))
  Rdat.append(tR / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Sdat)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
outfile = open('results/bilin.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
print >> outfile, "# diff err rel err gauge err SU(N) err U(1) err"

dat = np.array(Odat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Rdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Gdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Sdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> outfile, "%.6g %.4g" % (ave, err),

dat = np.array(Udat)
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
