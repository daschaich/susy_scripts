#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA susceptibility averages
# with blocked standard errors

# The C measurement now does all the A4* lattice conversion for us
# We just need to read, construct
#   \int dr r^{n + 3} \cO(x_0 + r) \cO(x_0) = \int dr r^{n + 3} C(r)
# for 0 <= n <= 3, and average

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

# Cycle through first output file to associate scalar distances
# with the corresponding label in the output
r = []          # List of r sorted by tag
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nt'):
    Nt = int((line.split())[1])
  # Format: CORR_K label r dat
  elif line.startswith('CORR_K '):
    r.append(float((line.split())[2]))
  elif line.startswith('CORR_S '):
    break             # Don't go through whole file yet

# Number of powers to consider: 0 <= n <= 4
Nsus = 4
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked susceptibilities for each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
Kdat = [[] for x in range(Nsus)]
Sdat = [[] for x in range(Nsus)]

# Make time series plot of operators themselves at the same time
timeseries = open('data/konishi.csv', 'w')
print >> timeseries, "MDTU,Konishi,SUGRA"

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [0 for x in range(Nsus)]
tS = [0 for x in range(Nsus)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for i in range(Nsus):
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
  OK = 0.0
  OS = 0.0
  check = -1
  for line in open(toOpen[0]):
    # Format: CORR_{K, S} label r dat
    if line.startswith('CORR'):
      temp = line.split()
      r = float(temp[2])
      dat = float(temp[3])
      if line.startswith('CORR_K 0 '):
        count += 1          # Only increment once per measurement!

      if r < 1.5:           # !!! Drop first three points
        continue

      if line.startswith('CORR_K '):
        for n in range(Nsus):
          tK[n] += np.power(r, n + 3) * dat
      elif line.startswith('CORR_S '):
        for n in range(Nsus):
          tS[n] += np.power(r, n + 3) * dat
      else:
        print "ERROR: tag ", temp[0], "not recognized"
        sys.exit(1)
    elif line.startswith('KONISHI '):
      OK += float((line.split())[2])
    elif line.startswith('SUGRA '):
      OS += float((line.split())[2])
    elif line.startswith('RUNNING COMPLETED'):
      # Average Konishi and SUGRA over all time slices
      OK /= Nt
      OS /= Nt
      print >> timeseries, "%d,%.6g,%.6g" % (MDTU, OK, OS)

      if check == 1:    # Check that we have one measurement per file
        print infile, "reports two measurements"
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)
timeseries.close()

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for i in range(Nsus):
    Kdat[i].append(tK[i] / float(count))
    Sdat[i].append(tS[i] / float(count))
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
outfile = open('results/suscept.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# n Konishi err SUGRA err"

for i in range(Nsus):
  print >> outfile, "%d" % i,

  dat = np.array(Kdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.6g %.4g" % (ave, err),

  dat = np.array(Sdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> outfile, "%.6g %.4g" % (ave, err)
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
