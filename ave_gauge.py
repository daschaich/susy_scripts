#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked gauge observable averages with blocked standard error
# For those ensembles I didn't generate
# (and so aren't in the format expected by average_SYM.py)

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
# F = Tr[Udag U] / N ('flink')
# P = plaquette, L = Polyakov loop, B = bosonic action
Fdat = []
Pdat = []
Ldat = []
Bdat = []

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tF = 0.0
tP = 0.0
tL = 0.0
tB = 0.0
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    Fdat.append(tF / float(count))
    Pdat.append(tP / float(count))
    Ldat.append(tL / float(count))
    Bdat.append(tB / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tF = 0.0
    tP = 0.0
    tL = 0.0
    tB = 0.0

  # Running average
  filename = 'Out/corr.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = 1
  for line in open(toOpen[0]):
    if line.startswith('FLINK '):
      # Format: FLINK link[0] link[1] link[2] link[3] link[4] ave
      count += 1    # Only increment once
      temp = line.split()
      tF += float(temp[6])
    elif line.startswith('GMES '):
      # Format: GMES Re(poly) Im(poly) junk ss_plaq st_plaq SB
      temp = line.split()
      tL += np.sqrt((float(temp[1]))**2 + (float(temp[2]))**2)
      tP += (float(temp[4]) + float(temp[5])) / 2.0
      tB += float(temp[6])
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  Fdat.append(tF / float(count))
  Pdat.append(tP / float(count))
  Ldat.append(tL / float(count))
  Bdat.append(tB / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Fdat)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Ffile = open('results/Flink.dat', 'w')
Pfile = open('results/plaq.dat', 'w')
Lfile = open('results/poly_mod.dat', 'w')
Bfile = open('results/SB.dat', 'w')
for outfile in [Ffile, Pfile, Lfile, Bfile]:
  print >> outfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
  print >> outfile, "# diff err rel err gauge err bilin err"

dat = np.array(Fdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> Ffile, "%.6g %.4g" % (ave, err),

dat = np.array(Pdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> Pfile, "%.6g %.4g" % (ave, err),

dat = np.array(Ldat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> Lfile, "%.6g %.4g" % (ave, err),

dat = np.array(Bdat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
print >> Bfile, "%.6g %.4g" % (ave, err),

for outfile in [Ffile, Pfile, Lfile, Bfile]:
  print >> outfile, "# %d" % Nblocks
  # More detailed block information
  #for i in range(Nblocks):
  #  print >> outfile, \
  #        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
  #        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
  outfile.close()
# ------------------------------------------------------------------
