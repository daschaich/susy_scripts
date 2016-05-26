#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA vacuum expectation values
# with blocked standard errors

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

numK = 3
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

# Cycle through first output file to find out the number of timeslices
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  # Format: CORR_K label r dat
  if line.startswith('nt '):
    Nt = float((line.split())[1])
    break             # Don't go through whole file yet
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each operator
# K = Konishi, S = SUGRA (averaged over all independent components)
# Check measurements with volume average subtracted, which should vanish
Kdat = [[] for i in range(2 * numK)]
Sdat = [[] for i in range(2 * numK)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = np.zeros(2 * numK, dtype = np.float)
tS = np.zeros_like(tK)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for i in range(2 * numK):
      Kdat[i].append(tK[i] / float(Nt * count))
      Sdat[i].append(tS[i] / float(Nt * count))
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
    # Format: TAG t op dat
    if line.startswith('KONISHI '):
      temp = line.split()
      op = int(temp[2])
      tK[op] += float(temp[3])
      tK[numK + op] += float(temp[4])
      if line.startswith('KONISHI 0 0 '):
        count += 1          # Only increment once per measurement!

    elif line.startswith('SUGRA '):
      temp = line.split()
      op = int(temp[2])
      tS[op] += float(temp[3])
      tS[numK + op] += float(temp[4])
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
  for i in range(2 * numK):
    Kdat[i].append(tK[i] / float(Nt * count))
    Sdat[i].append(tS[i] / float(Nt * count))
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
outfile = open('results/vevK.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# Konishi err"

for i in range(2 * numK):
  dat = np.array(Kdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.16g %.4g" % (ave, err)

print >> outfile, "# SUGRA err"
for i in range(2 * numK):
  dat = np.array(Sdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.16g %.4g" % (ave, err)
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
