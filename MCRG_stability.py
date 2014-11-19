#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Compute MCRG stability matrix using blocked jackknife procedure
# Compute xi^4 by matching blocked plaquette

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us which files to analyze
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

# Determine maximum blocking level from first output file
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
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
# Construct arrays of blocked and unblocked measurements
# for the Konishi (K) and SUGRA (S) operators
# We also need to accumulate two products for each operator
# We also need (blmax - 1) 2x2 and 1x1 Wilson loops to determine xi^4
Kdat = [[] for i in range(blmax + 1)]
AKdat = [[] for i in range(blmax)]
BKdat = [[] for i in range(blmax)]
Sdat = [[] for i in range(blmax + 1)]
ASdat = [[] for i in range(blmax)]
BSdat = [[] for i in range(blmax)]
L2dat = [[] for i in range(blmax)]
L1dat = [[] for i in range(blmax)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = np.zeros(blmax + 1, dtype = np.float)
tAK = np.zeros(blmax, dtype = np.float)
tBK = np.zeros(blmax, dtype = np.float)
tS = np.zeros(blmax + 1, dtype = np.float)
tAS = np.zeros(blmax, dtype = np.float)
tBS = np.zeros(blmax, dtype = np.float)
tL2 = np.zeros(blmax, dtype = np.float)
tL1 = np.zeros(blmax, dtype = np.float)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    if count == 0:
      print "ERROR: no data to average after file %s:" % toOpen
      sys.exit(1)
    for bl in range(blmax):
      Kdat[bl].append(tK[bl] / float(count))
      AKdat[bl].append(tAK[bl] / float(count))
      BKdat[bl].append(tBK[bl] / float(count))
      Sdat[bl].append(tS[bl] / float(count))
      ASdat[bl].append(tAS[bl] / float(count))
      BSdat[bl].append(tBS[bl] / float(count))
      L2dat[bl].append(tL2[bl] / float(20. * count))
      L1dat[bl].append(tL2[bl] / float(20. * count))
      tK[bl] = 0.0
      tAK[bl] = 0.0
      tBK[bl] = 0.0
      tS[bl] = 0.0
      tAS[bl] = 0.0
      tBS[bl] = 0.0
      tL2[bl] = 0.0
      tL1[bl] = 0.0
    Kdat[blmax].append(tK[blmax] / float(count))
    Sdat[blmax].append(tS[blmax] / float(count))
    tK[blmax] = 0.0
    tS[blmax] = 0.0
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
  for line in open(toOpen[0]):
    # Format: OK bl blocked_konishi
    if line.startswith('OK '):
      temp = line.split()
      bl = int(temp[1])
      dat = float(temp[2])
      tK[bl] += dat

      # Accumulate products as well as operator -- note shifted index
      # B(n) is easy, A(n) requires saving smaller blocking level
      if bl > 0:
        tBK[bl - 1] = dat * dat
        if bl < blmax:
          tAK[bl - 1] = dat * savK
          savK = dat
      if bl == 0:
        savK = dat
        count += 1          # Only tick counter once per measurement

    # Format: OS bl blocked_SUGRA
    elif line.startswith('OS '):
      temp = line.split()
      bl = int(temp[1])
      dat = float(temp[2])
      tS[bl] += dat

      # Accumulate products as well as operator -- note shifted index
      # B(n) is easy, A(n) requires saving smaller blocking level
      if bl > 0:
        tBS[bl - 1] = dat * dat
        if bl < blmax:
          tAS[bl - 1] = dat * savS
          savS = dat
      if bl == 0:
        savS = dat

    # Format: BRSYMM bl normal [dir] inverted [dir] usual mod
    elif line.startswith('BRSYMM '):
      temp = line.split()
      bl = int(temp[1])
      norm = int(temp[2])
      inv = int(temp[4])
      # 2x2 loops on non-maximally blocked lattices
      if bl < blmax and norm == 2 and inv == 2:
        tL2[bl] += float(temp[6])
      # 1x1 loops on blocked lattices -- note shifted index
      elif bl > 0 and norm == 1 and inv == 1:
        tL1[bl - 1] += float(temp[6])

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  if count == 0:
    print "ERROR: no data to average after file %s:" % toOpen
    sys.exit(1)
  for bl in range(blmax):
    Kdat[bl].append(tK[bl] / float(count))
    Sdat[bl].append(tS[bl] / float(count))
    L2dat[bl].append(tL2[bl] / float(20. * count))
    L1dat[bl].append(tL1[bl] / float(20. * count))
  Kdat[blmax].append(tK[blmax] / float(count))
  Sdat[blmax].append(tS[blmax] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination,
# and analyze them to determine xi^4(n), A(n), B(n) and T(n),
# eventually obtaining jackknife estimates for Delta = 4 - y
# for both Konishi and SUGRA
# Require multiple blocks instead of worrying about error propagation
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

Ktot = np.empty_like(tK)
AKtot = np.empty_like(tL2)
BKtot = np.empty_like(tL2)
Stot = np.empty_like(tK)
AStot = np.empty_like(tL2)
BStot = np.empty_like(tL2)
L2tot = np.empty_like(tL2)
L1tot = np.empty_like(tL1)
for bl in range(blmax):
  Ktot[bl] = sum(Kdat[bl])
  AKtot[bl] = sum(AKdat[bl])
  BKtot[bl] = sum(BKdat[bl])
  Stot[bl] = sum(Sdat[bl])
  AStot[bl] = sum(ASdat[bl])
  BStot[bl] = sum(BSdat[bl])
  L2tot[bl] = sum(L2dat[bl])
  L1tot[bl] = sum(L1dat[bl])
Ktot[blmax] = sum(Kdat[blmax])
Stot[blmax] = sum(Sdat[blmax])

# All jackknife results
# Only care about Delta, but might be worthwhile to monitor xi^4
jkXi = np.zeros((blmax, Nblocks), dtype = np.float)
jkDeltaK = np.zeros((blmax, Nblocks), dtype = np.float)
jkDeltaS = np.zeros_like(jkDeltaK)

for i in range(Nblocks):  # Jackknife samples
  # Compute xi^4 for each blocking level
  for bl in range(blmax):
    L1 = L1tot[bl] - L1dat[bl][i]
    L2 = L2tot[bl] - L2dat[bl][i]
    print bl, L1, L2
    jkXi[bl][i] = L1 / L2

    # It's a little awkward that the indices of the individual operators
    # are shifted relative to all the others
    temp = (Ktot[bl] - Kdat[bl][i]) / (Nblocks - 1.0)
    if bl > 0:
      temp *= jkXi[bl - 1][i]
    K = jkXi[bl][i] * (Ktot[bl + 1] - Kdat[bl + 1][i]) / (Nblocks - 1.0)
    AK = jkXi[bl][i] * (AKtot[bl] - AKdat[bl][i]) / (Nblocks - 1.0)
    AK -= K * temp    # Subtract product of estimates
    BK = jkXi[bl][i] * (BKtot[bl] - BKdat[bl][i]) / (Nblocks - 1.0)
    BK -= K * K       # Subtract product of estimates
    TK = AK / BK
    jkDeltaK[bl][i] = np.log(TK) / np.log(2.0)

    temp = (Stot[bl] - Sdat[bl][i]) / (Nblocks - 1.0)
    if bl > 0:
      temp *= jkXi[bl - 1][i]
    S = jkXi[bl][i] * (Stot[bl + 1] - Sdat[bl + 1][i]) / (Nblocks - 1.0)
    AS = jkXi[bl][i] * (AStot[bl] - ASdat[bl][i]) / (Nblocks - 1.0)
    AS -= S * temp    # Subtract product of estimates
    BS = jkXi[bl][i] * (BStot[bl] - BSdat[bl][i]) / (Nblocks - 1.0)
    BS -= S * S       # Subtract product of estimates
    TS = AS / BS
    jkDeltaS[bl][i] = np.log(TS) / np.log(2.0)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)
outfile = open('results/MCRG.dat', 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)

for bl in range(blmax):
  ave = np.average(jkXi[bl])
  err = (Nblocks - 1.0) * np.sum((jkXi[bl] - ave)**2) / float(Nblocks)
  print >> outfile, "xi^4 %d %.6g %.4g" % (bl, ave, err)

  ave = np.average(jkDeltaK[bl])
  err = (Nblocks - 1.0) * np.sum((jkDeltaK[bl] - ave)**2) / float(Nblocks)
  print >> outfile, "DeltaK %d %.6g %.4g" % (bl, ave, err)

  ave = np.average(jkDeltaS[bl])
  err = (Nblocks - 1.0) * np.sum((jkDeltaS[bl] - ave)**2) / float(Nblocks)
  print >> outfile, "DeltaS %d %.6g %.4g" % (bl, ave, err)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
