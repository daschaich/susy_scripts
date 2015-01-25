#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Blocked Konishi and SUGRA correlators C(t) and Delta C / Delta t
# projected to zero spatial momentum

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
# To only print 0 through Nt / 2 (a total of Npts points)
Npts = Nt / 2 + 1  # Assume Nt is even
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA
# DK and DS are corresponding finite differences
Kdat  = [[] for x in range(Npts)]
DKdat = [[] for x in range(Npts - 1)]
Sdat  = [[] for x in range(Npts)]
DSdat = [[] for x in range(Npts - 1)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK  = [0 for x in range(Npts)]
tDK = [0 for x in range(Npts - 1)]
tS  = [0 for x in range(Npts)]
tDS = [0 for x in range(Npts - 1)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for t in range(Npts):
      Kdat[t].append(tK[t] / float(count))
      Sdat[t].append(tS[t] / float(count))
    for t in range(1, Npts):
      DKdat[t - 1].append(tDK[t - 1] / float(count))
      DSdat[t - 1].append(tDS[t - 1] / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK  = [0 for x in range(Npts)]
    tDK = [0 for x in range(Npts - 1)]
    tS  = [0 for x in range(Npts)]
    tDS = [0 for x in range(Npts - 1)]

  # Running averages
  prev_time = float('NaN')    # To make it obvious if this isn't overwritten
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  for line in open(toOpen[0]):
    # Format: KONISHI t dat
    if line.startswith('KONISHI '):
      temp = line.split()
      t = int(temp[1])
      dat = float(temp[2])
      tK[t] += dat
      if t == 0:
        count += 1    # Only increment once per measurement!
      elif t > 0:
        tDK[t - 1] += dat - prev_time
      prev_time = dat
    # We go through all Konishi data before reaching first SUGRA line
    # Format: SUGRA t dat
    elif line.startswith('SUGRA '):
      temp = line.split()
      t = int(temp[1])
      dat = float(temp[2])
      tS[t] += dat
      if t > 0:
        tDS[t - 1] += dat - prev_time
      prev_time = dat
    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print infile, "reports two measurements"
        print >> ERRFILE, infile, "reports two measurements"
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for t in range(Npts):
    Kdat[t].append(tK[t] / float(count))
    Sdat[t].append(tS[t] / float(count))
  for t in range(1, Npts):
    DKdat[t - 1].append(tDK[t - 1] / float(count))
    DSdat[t - 1].append(tDS[t - 1] / float(count))
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
Kfile = open('results/konishi_t.dat', 'w')
print >> Kfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Sfile = open('results/sugra_t.dat', 'w')
print >> Sfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
for t in range(Npts):
  # Konishi
  dat = np.array(Kdat[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> Kfile, "%d %.6g %.4g" % (t, ave, err)

  # SUGRA
  dat = np.array(Sdat[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> Sfile, "%d %.6g %.4g" % (t, ave, err)
Kfile.close()
Sfile.close()

DKfile = open('results/konishi_diff.dat', 'w')
print >> DKfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
DSfile = open('results/sugra_diff.dat', 'w')
print >> DSfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
for t in range(1, Npts):
  # Konshi finite difference in time
  dat = np.array(DKdat[t - 1])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> DKfile, "%.2g %.6g %.4g" % (t - 0.5, ave, err)

  # Konshi finite difference in time
  dat = np.array(DSdat[t - 1])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> DSfile, "%.2g %.6g %.4g" % (t - 0.5, ave, err)
DKfile.close()
DSfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
