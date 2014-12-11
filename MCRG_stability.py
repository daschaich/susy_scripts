#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Compute MCRG stability matrix using blocked jackknife procedure
# Consider up to three operators per measurement
# Limitation: The order Tr(BB) Tr(BBBB) Tr(BB)Tr(BB) is fixed

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us how many operators to consider per measurement
# Fourth argument tells us how many directories to analyze,
# which are listed as any remaining arguments
# Each of these directories must contain mcrg files,
# and optionally xi^4 files to load
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <# ops> <# dirs> <dirs>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
ops_per_dir = int(sys.argv[3])
num_dirs = int(sys.argv[4])
num_ops = ops_per_dir * num_dirs
dirs = []
for i in range(num_dirs):
  dirs.append(str(sys.argv[5 + i]))
runtime = -time.time()

# Quick sanity check
if ops_per_dir < 1 or ops_per_dir > 3:
  print "ERROR: Only have 1, 2 or 3 operators per measurement"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Use first directory to construct list of analyzed configurations
cfgs = []
files = dirs[0] + '/mcrg.*'
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

# Determine maximum blocking level from first output file
firstfile = dirs[0] + '/mcrg.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)

blmax = -1
for line in open(firstfile):
  if line.startswith('Blocking '):
    blmax = int((line.split())[1])

if blmax < 0:
  print "ERROR:", firstfile, "doesn't mention blocking"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Either load xi array from files (ignoring uncertainties) or set it to unity
xi = np.zeros((num_ops, blmax + 1), dtype = np.float)
for i in range(num_dirs):
  xi_file = dirs[i] + '/xi.dat'
  if os.path.isfile(xi_file):
    print "Reading xi^4 from", xi_file
    for line in open(xi_file):
      if line.startswith('# '):
        continue
      else:
        temp = line.split()
        bl = int(temp[0])
        for j in range(ops_per_dir):
          xi[ops_per_dir * i + j][bl] = float(temp[1])
  else:
    print "Setting xi^4 to unity for directory", dirs[i]
    for bl in range(blmax + 1):
      for j in range(ops_per_dir):
        xi[ops_per_dir * i + j][bl] = 1.0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements
# for the Konishi (K) and SUGRA (S) operators
# We also need to accumulate two arrays of products for each operator
# These list definitions are a little awkward: the last index comes first
Kdat = [[[] for i in range(blmax + 1)] for j in range(num_ops)]
AKdat = [[[[] for i in range(blmax)] for j in range(num_ops)] for k in range(num_ops)]
BKdat = [[[[] for i in range(blmax)] for j in range(num_ops)] for k in range(num_ops)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = np.zeros((num_ops, blmax + 1), dtype = np.float)
tAK = np.zeros((num_ops, num_ops, blmax), dtype = np.float)
tBK = np.zeros((num_ops, num_ops, blmax), dtype = np.float)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    if count == 0:
      print "ERROR: no data to average after file %s:" % toOpen
      sys.exit(1)
    for i in range(num_ops):
      Kdat[i][blmax].append(tK[i][blmax] / float(count))
      tK[i][blmax] = 0.0
      for bl in range(blmax):
        Kdat[i][bl].append(tK[i][bl] / float(count))
        tK[i][bl] = 0.0
        for j in range(num_ops):
          AKdat[i][j][bl].append(tAK[i][j][bl] / float(count))
          BKdat[i][j][bl].append(tBK[i][j][bl] / float(count))
          tAK[i][j][bl] = 0.0
          tBK[i][j][bl] = 0.0
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)

  # Running averages require data from every directory
  dat = np.zeros((num_ops, blmax + 1), dtype = np.float)
  for i in range(num_dirs):
    filename = dirs[i] + '/mcrg.' + str(MDTU)
    toOpen = glob.glob(filename)
    if len(toOpen) > 1:
      print "ERROR: multiple files named %s:" % filename,
      print toOpen
    check = -1
    for line in open(toOpen[0]):
      # Format: OK bl Tr(BB) Tr(BBB) Tr(BBBB) Tr(BB)Tr(BB)
      if line.startswith('OK '):
        temp = line.split()
        bl = int(temp[1])
        dat[ops_per_dir * i][bl] = float(temp[2])
        if ops_per_dir > 1:
          # Tr(BBB) in temp[3] always vanishes
          dat[ops_per_dir * i + 1][bl] = float(temp[4])
        if ops_per_dir > 2:
          dat[ops_per_dir * i + 2][bl] = float(temp[5])
      elif line.startswith('RUNNING COMPLETED'):
        check = 1
    if check == -1:
      print toOpen[0], "did not complete"
      sys.exit(1)

  # Accumulate operator and products A(n), B(n) -- note shifted index
  count += 1                        # Only tick once per measurement
  for i in range(num_ops):
    tK[i][blmax] += dat[i][blmax]
    for bl in range(blmax):
      tK[i][bl] += dat[i][bl]
      for j in range(num_ops):
        tAK[i][j][bl] += dat[i][bl + 1] * dat[j][bl]
        tBK[i][j][bl] += dat[i][bl + 1] * dat[j][bl + 1]

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  if count == 0:
    print "ERROR: no data to average after file %s:" % toOpen
    sys.exit(1)
  for i in range(num_ops):
    Kdat[i][blmax].append(tK[i][blmax] / float(count))
    for bl in range(blmax):
      Kdat[i][bl].append(tK[i][bl] / float(count))
      for j in range(num_ops):
        AKdat[i][j][bl].append(tAK[i][j][bl] / float(count))
        BKdat[i][j][bl].append(tBK[i][j][bl] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination,
# and analyze them to determine xi^4(n), A(n), B(n) and T(n),
# eventually obtaining jackknife estimates for Delta = 4 - y
# for both Konishi and SUGRA
# Require multiple blocks instead of worrying about error propagation
if Nblocks == 1:
  print "ERROR: need multiple blocks to analyze"
  sys.exit(1)

Ktot = np.empty_like(tK)
AKtot = np.empty_like(tAK)
BKtot = np.empty_like(tBK)
for i in range(num_ops):
  Ktot[i][blmax] = sum(Kdat[i][blmax])
  for bl in range(blmax):
    Ktot[i][bl] = sum(Kdat[i][bl])
    for j in range(num_ops):
      AKtot[i][j][bl] = sum(AKdat[i][j][bl])
      BKtot[i][j][bl] = sum(BKdat[i][j][bl])

# All jackknife results -- only care about Delta = 4 - y
jkDeltaK = np.zeros((blmax, Nblocks), dtype = np.float)
for n in range(Nblocks):  # Jackknife samples
  # It's a little awkward that the indices of the individual operators
  # are shifted relative to all the others
  K = np.zeros(num_ops, dtype = np.float)
  temp = np.zeros_like(K)
  AK = np.zeros((num_ops, num_ops), dtype = np.float)
  BK = np.zeros_like(AK)
  for bl in range(blmax):
    for i in range(num_ops):
      temp[i] = (Ktot[i][bl] - Kdat[i][bl][n]) / (Nblocks - 1.0)
      temp[i] *= xi[i][bl]
      K[i] = (Ktot[i][bl + 1] - Kdat[i][bl + 1][n]) / (Nblocks - 1.0)
      K[i] *= xi[i][bl + 1]
    for i in range(num_ops):
      for j in range(num_ops):
        AK[i][j] = (AKtot[i][j][bl] - AKdat[i][j][bl][n]) / (Nblocks - 1.0)
        AK[i][j] *= xi[i][bl + 1] * xi[j][bl]
        AK[i][j] -= K[i] * temp[j]

        BK[i][j] = (BKtot[i][j][bl] - BKdat[i][j][bl][n]) / (Nblocks - 1.0)
        BK[i][j] *= xi[i][bl + 1] * xi[j][bl + 1]
        BK[i][j] -= K[i] * K[j]

    BKinv = np.linalg.inv(BK)
    TK = np.dot(BKinv, AK)
    eig, vecs = np.linalg.eig(TK)
    la_max = np.amax(np.absolute(eig))    # Absolute values...
    jkDeltaK[bl][n] = 4.0 - np.log(la_max) / np.log(2.0)

    # Optional monitoring for initial debugging / sanity checks
#    if n == Nblocks - 1:
#      print "AK:", AK
#      print "BK:", BK
#      print "BKinv:", BKinv
#      print "TK:", TK
#      print eig, "-->", la_max
#      print ""
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)
print "%d operators in the stability matrix" % num_ops
outfile = open('results/MCRG.dat', 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)
print >> outfile, "# %d operators in the stability matrix" % num_ops

for bl in range(blmax):
  ave = np.average(jkDeltaK[bl])
  err = (Nblocks - 1.0) * np.sum((jkDeltaK[bl] - ave)**2) / float(Nblocks)
  print >> outfile, "DeltaK %d %.6g %.4g" % (bl + 1, ave, err)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
