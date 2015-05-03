#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA operator averages
# with blocked standard errors

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us which files to analyze -- needs dir path
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <dir/tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct list of which configurations have been analyzed
cfgs = []
files = tag + '.*'
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

# Determine largest Wilson loops from first output file
firstfile = tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)

# Extract volume and number of blocking levels
vol = -1;
blmax = -1
for line in open(firstfile):
  if line.startswith('Blocking '):
    blmax = int((line.split())[1])
  elif line.startswith('nx '):
    vol = float((line.split())[1])
  elif line.startswith('ny ') or line.startswith('nz '):
    vol *= float((line.split())[1])
  elif line.startswith('nt '):
    vol *= float((line.split())[1])

if vol < 0:
  print "ERROR:", firstfile, "doesn't define lattice volume"
  sys.exit(1)
if blmax < 0:
  print "ERROR:", firstfile, "doesn't mention blocking"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each observable
# for the Konishi (K) and SUGRA (S) operators
Kdat = [[] for i in range(blmax + 1)]
Sdat = [[] for i in range(blmax + 1)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = np.zeros(blmax + 1, dtype = np.float)
tS = np.zeros(blmax + 1, dtype = np.float)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    if count == 0:
      print "ERROR: no data to average after file %s:" % toOpen
      sys.exit(1)
    for bl in range(blmax + 1):
      Kdat[bl].append(tK[bl] / float(count * vol / 16.0**bl))
      Sdat[bl].append(tS[bl] / float(count * vol / 16.0**bl))
      tK[bl] = 0.0
      tS[bl] = 0.0
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)

  # Running averages
  filename = tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  for line in open(toOpen[0]):
    # Format: OK smearing bl blocked_konishi
    if line.startswith('OK '):
      temp = line.split()
      bl = int(temp[2])
      tK[bl] += float(temp[3])
      if bl == 0:
        count += 1          # Only tick counter once per measurement
    # Format: OS bl blocked_SUGRA
    elif line.startswith('OS '):
      temp = line.split()
      bl = int(temp[2])
      tS[bl] += float(temp[3])
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  if count == 0:
    print "ERROR: no data to average after file %s:" % toOpen
    sys.exit(1)
  for bl in range(blmax + 1):
    Kdat[bl].append(tK[bl] / float(count * vol / 16.0**bl))
    Sdat[bl].append(tS[bl] / float(count * vol / 16.0**bl))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
outfile = open('results/ops.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
print >> outfile, "# bl Konishi err SUGRA err"

for i in range(blmax + 1):
  print >> outfile, i,

  dat = np.array(Kdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> outfile, "%.6g %.4g" % (ave, err),

  dat = np.array(Sdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> outfile, "%.6g %.4g" % (ave, err)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
