#!/usr/bin/python
import os
import sys
import glob
import time
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
runtime = -time.time()

# Convenience constants for sanity check
TOL = 1.0e-6
invSq2  = 1.0 / np.sqrt(2)
invSq6  = 1.0 / np.sqrt(6)
invSq12 = 1.0 / np.sqrt(12)
invSq20 = 1.0 / np.sqrt(20)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Map (x, y, z, t) to r on L^3 x Nt A4* lattice
# Check all possible periodic shifts to find true r
def A4map(x_in, y_in, z_in, t_in, L, Nt):
  r = np.sqrt((x_in**2 + y_in**2 + z_in**2 + t_in**2) * 0.8 \
              - (x_in * (y_in + z_in + t_in) \
                 + y_in * (z_in + t_in) + z_in * t_in) * 0.4)
  for x in [x_in + L, x_in, x_in - L]:
    for y in [y_in + L, y_in, y_in - L]:
      for z in [z_in + L, z_in, z_in - L]:
        for t in [t_in + Nt, t_in, t_in - Nt]:
          test = np.sqrt((x**2 + y**2 + z**2 + t**2) * 0.8 \
                         - (x * (y + z + t) \
                            + y * (z + t) + z * t) * 0.4)

          # Sanity check -- can be commented out for more speed
#          x_a4 = (x - y) * invSq2
#          y_a4 = (x + y - 2.0 * z) * invSq6
#          z_a4 = (x + y + z - 3.0 * t) * invSq12
#          t_a4 = (x + y + z + t) * invSq20
#          check = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2 + t_a4**2)
#          if np.fabs(test - check) > TOL:
#            print "ERROR: %.6g isn't %.6g for (%d, %d, %d, %d)" \
#                  % (check, test, x, y, z, t)
#            sys.exit(1)

          # True r is the smallest
          # Try to avoid negative roundoff
          # so that we can see when periodic shifts are really necessary
          if test - r < -TOL:
#            print "|(%d, %d, %d, %d)| = %.6g" % (x_in, y_in, z_in, t_in, r),
#            print "--> |(%d, %d, %d, %d)| = %.6g" % (x, y, z, t, test)
            r = test
#  print "%d %d %d %d --> %.6g" % (x_in, y_in, z_in, t_in, r)
  return r
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
# To simplify things later, record these scalar distances in a lookup table
# Also, decide where to cut r based on MAX_X
r = []          # List of two-component lists: first value, then count
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('d_correlator_r: MAX_T = '):
    MAX_X = int((line.split())[6])
    num_y = 2 * MAX_X + 1
    lookup = np.empty((MAX_X + 1, num_y, num_y, num_y), dtype = np.float)

    # Decide where to cut r
    MAX_r = 100.0 * MAX_X   # To be overwritten
    for y in range(MAX_X + 1):
      for z in range(MAX_X + 1):
        for t in range(MAX_X + 1):
          test = A4map(MAX_X + 1, y, z, t, L, Nt)
          if test < MAX_r:
            MAX_r = test
#    print MAX_r

  # Format: CORR_K x y z t dat
  elif line.startswith('CORR_K '):
    temp = line.split()
    x = int(temp[1])
    y = int(temp[2])
    z = int(temp[3])
    t = int(temp[4])
    if t > MAX_X:         # Drop large separations
      continue
    this_r = A4map(x, y, z, t, L, Nt)

    # Save this_r in lookup table so we can avoid calling A4map in the future
    # Note that we require t no larger than MAX_X
    lookup[x][y + MAX_X][z + MAX_X][t + MAX_X] = this_r

    # Accumulate multiplicities if this_r < MAX_r
    # Try to avoid roundoff issues in comparison
    if this_r - MAX_r > -TOL:
      continue
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
#print lookup
#for j in range(Npts):
#  print r[j]
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
    # Format: CORR_[K,S] x y z t dat
    if line.startswith('CORR'):
      temp = line.split()
      x = int(temp[1])
      y = int(temp[2])
      z = int(temp[3])
      t = int(temp[4])
      dat = float(temp[5])
      if t > MAX_X:         # Drop large separations, as assumed by lookup
        continue
      if line.startswith('CORR_K 0 0 0 0 '):
        count += 1          # Only increment once per measurement!
      this_r = lookup[x][y + MAX_X][z + MAX_X][t + MAX_X]

      # Accumulate data if this_r < MAX_r
      # Try to avoid roundoff shenanigans
      if this_r - MAX_r > -TOL:
        continue
      done = -1
      for i in range(len(r)):
        if np.fabs(r[i][0] - this_r) < 1.0e-6:
          if line.startswith('CORR_K '):
            tK[i] += dat
          elif line.startswith('CORR_S '):
            tS[i] += dat
          else:
            print "ERROR: tag ", temp[0], "not recognized"
            sys.exit(1)
          done = 1
          break
      if done < 0:
        print "ERROR: displacement", this_r, "not found"
        sys.exit(1)
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
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
Kfile = open('results/konishi_r.dat', 'w')
print >> Kfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Sfile = open('results/sugra_r.dat', 'w')
print >> Sfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
for i in range(Npts):
  # Konishi
  dat = np.array(Kdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> Kfile, "%.4g %.6g %.4g" % (r[i][0], ave, err)

  # SUGRA
  dat = np.array(Sdat[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print >> Sfile, "%.4g %.6g %.4g" % (r[i][0], ave, err)
Kfile.close()
Sfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
