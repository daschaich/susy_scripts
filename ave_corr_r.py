#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA correlator averages
# with blocked standard errors

# The C measurement now does all the A4* lattice conversion for us
# We just need to read and average

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
r = []          # List of r sorted by tag
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
      r.append(float((line.split())[2]))
  elif line.startswith('CORR_S '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
# First subtracts ensemble average, second volume average
Kdat = [[[] for x in range(Npts)] for i in range(2)]
Sdat = [[[] for x in range(Npts)] for i in range(2)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [[[] for x in range(Npts)] for i in range(2)]
tS = [[[] for x in range(Npts)] for i in range(2)]
for sub in range(2):
  for i in range(Npts):
    tK[sub][i] = 0.0
    tS[sub][i] = 0.0

for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for sub in range(2):
      for i in range(Npts):
        Kdat[sub][i].append(tK[sub][i] / float(count))
        Sdat[sub][i].append(tS[sub][i] / float(count))
        tK[sub][i] = 0.0
        tS[sub][i] = 0.0
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
      temp = line.split()
      label = int(temp[1])
      tag1 = int(temp[3])
      tag2 = int(temp[4])
      dat_vev = float(temp[5])
      dat_vol = float(temp[6])
      if line.startswith('CORR_K 0 0 0 0 '):
        count += 1          # Only increment once per measurement!

      if line.startswith('CORR_K ') and tag1 == 0 and tag2 == 0:
        tK[0][label] += dat_vev
        tK[1][label] += dat_vol
      elif line.startswith('CORR_S ') and tag1 == 0 and tag2 == 0:
        tS[0][label] += dat_vev
        tS[1][label] += dat_vol
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
  for sub in range(2):
    for i in range(Npts):
      Kdat[sub][i].append(tK[sub][i] / float(count))
      Sdat[sub][i].append(tS[sub][i] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

columns = [('r', float), ('ave', float), ('err', float)]
Kout_vev = np.zeros(Npts, dtype = columns)
Kout_vol = np.zeros_like(Kout_vev)
Sout_vev = np.zeros_like(Kout_vev)
Sout_vol = np.zeros_like(Kout_vev)
for i in range(Npts):
  # Konishi with ensemble average subtracted
  dat = np.array(Kdat[0][i])
  Kout_vev[i][0] = r[i]
  Kout_vev[i][1] = np.mean(dat, dtype = np.float64)
  Kout_vev[i][2] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)

  # Konishi with volume average subtracted
  dat = np.array(Kdat[1][i])
  Kout_vol[i][0] = r[i]
  Kout_vol[i][1] = np.mean(dat, dtype = np.float64)
  Kout_vol[i][2] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)

  # SUGRA with ensemble average subtracted
  dat = np.array(Sdat[0][i])
  Sout_vev[i][0] = r[i]
  Sout_vev[i][1] = np.mean(dat, dtype = np.float64)
  Sout_vev[i][2] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)

  # SUGRA with volume average subtracted
  dat = np.array(Sdat[1][i])
  Sout_vol[i][0] = r[i]
  Sout_vol[i][1] = np.mean(dat, dtype = np.float64)
  Sout_vol[i][2] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)

# Sort output by r (column zero) and print
Kout_vev = np.sort(Kout_vev, order='r')
Kout_vol = np.sort(Kout_vol, order='r')
Sout_vev = np.sort(Sout_vev, order='r')
Sout_vol = np.sort(Sout_vol, order='r')

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Kfile = open('results/konishi_r.dat', 'w')
print >> Kfile, "# Averaging with %d blocks of length %d MDTU" \
                % (Nblocks, block_size)
print >> Kfile, "# r    vev       err       vol       err"
Sfile = open('results/sugra_r.dat', 'w')
print >> Sfile, "# Averaging with %d blocks of length %d MDTU" \
                % (Nblocks, block_size)
print >> Sfile, "# r    vev       err       vol       err"
for i in range(Npts):
  print >> Kfile, "%.4g %.6g %.4g" \
                  % (Kout_vev[i][0], Kout_vev[i][1], Kout_vev[i][2]),
  print >> Kfile, "%.6g %.4g" % (Kout_vol[i][1], Kout_vol[i][2])
  print >> Sfile, "%.4g %.6g %.4g" \
                  % (Sout_vev[i][0], Sout_vev[i][1], Sout_vev[i][2]),
  print >> Sfile, "%.6g %.4g" % (Sout_vol[i][1], Sout_vol[i][2])
Kfile.close()
Sfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
