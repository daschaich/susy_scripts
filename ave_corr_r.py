#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA correlator averages
# with blocked standard errors

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

# Convenience constants for sanity check
invSq2  = 1.0 / np.sqrt(2)
invSq6  = 1.0 / np.sqrt(6)
invSq12 = 1.0 / np.sqrt(12)
invSq20 = 1.0 / np.sqrt(20)
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

# Cycle through first output file to construct list of scalar distances
# and determine their multiplicity
r = []          # List of two-component lists: first value, then count
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  # Format: CORR_K x y z t dat
  if line.startswith('CORR_K '):
    temp = line.split()
    x = float(temp[1])
    y = float(temp[2])
    z = float(temp[3])
    t = float(temp[4])
    this_r = np.sqrt((x**2 + y**2 + z**2 + t**2) * 0.8 \
                     - 2.0 * (x * (y + z + t) + y * (z + t) + z * t) * 0.2)

    # Sanity check
    x_a4 = (x - y) * invSq2
    y_a4 = (x + y - 2.0 * z) * invSq6
    z_a4 = (x + y + z - 3.0 * t) * invSq12
    t_a4 = (x + y + z + t) * invSq20
    check_r = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2 + t_a4**2)
    if np.fabs(this_r - check_r) > 1.0e-6:
      print "ERROR:", this_r, "isn't", check_r
      sys.exit(1)

    # Inefficient, but whatever
    done = -1
    for i in range(len(r)):
      if np.fabs(r[i][0] - this_r) < 1.0e-6:
        r[i][1] += 1
        done = 1
        break
    if done < 0:
      r.append([this_r, 1])
Npts = len(r)

# Sort by magnitude (column zero), not count
r = sorted(r, key=lambda x: x[0])
#for i in range(Npts):
#  print r[i]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA (averaged over all 25 components)
Kdat = [[] for x in range(Npts)]
Sdat = [[] for x in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [0 for x in range(Npts)]
tS = [0 for x in range(Npts)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for i in range(Npts):
      Kdat[i].append(tK[i] / float(count * r[i][1]))
      Sdat[i].append(tS[i] / float(count * r[i][1]))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK = [0 for x in range(Npts)]
    tS = [0 for x in range(Npts)]

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: CORR_[K,S] x y z t dat
    if line.startswith('CORR'):
      temp = line.split()
      x = float(temp[1])
      y = float(temp[2])
      z = float(temp[3])
      t = float(temp[4])
      this_r = np.sqrt((x**2 + y**2 + z**2 + t**2) * 0.8 \
                       - 2.0 * (x * (y + z + t) + y * (z + t) + z * t) * 0.2)
      if this_r < 1.0e-6:   # (0, 0, 0, 0) only
        count += 1          # Only increment once per measurement!

      # Inefficient, but whatever
      done = -1
      for i in range(len(r)):
        if np.fabs(r[i][0] - this_r) < 1.0e-6:
          if line.startswith('CORR_K '):
            tK[i] += float(temp[5])
          elif line.startswith('CORR_S '):
            tS[i] += float(temp[5])
          else:
            print "ERROR: tag ", temp[0], "not recognized"
            sys.exit(1)
          done = 1
          break
      if done < 0:
        print "ERROR: displacement", this_r, "not found"
        sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for i in range(Npts):
    Kdat[i].append(tK[i] / float(count * r[i][1]))
    Sdat[i].append(tS[i] / float(count * r[i][1]))
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
Kfile = open('results/konishi.dat', 'w')
print >> Kfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Sfile = open('results/sugra.dat', 'w')
print >> Sfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
for i in range(Npts):
  # Konishi
  dat = np.array(Kdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> Kfile, "%.4g %.6g %.4g # %d" % (r[i][0], ave, err, Nblocks)

  # SUGRA
  dat = np.array(Sdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print >> Sfile, "%.4g %.6g %.4g # %d" % (r[i][0], ave, err, Nblocks)

# More detailed block information
#for i in range(Nblocks):
#  for outfile in [Kfile, Sfile]:
#    print >> outfile, \
#          "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#          % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
Kfile.close()
Sfile.close()
# ------------------------------------------------------------------
