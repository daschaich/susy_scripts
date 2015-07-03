#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Compute MCRG stability matrix using blocked jackknife procedure
# Consider only two operators per measurement
#   Tr(PP) is built from polar projections of the links -- give it xi^2
#   Tr(UU) is built from U.Ubar -- give it xi^4
# The order Tr(PP) Tr(UU) is fixed
# While lines around 145 can be adjusted below,
# also remember to adjust corresponding xi factors....

# Optionally load results/xi.dat, which I think should not depend on smearing
# since it is part of the RG blocking
# while the smearing is part of the observable definition

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us how many operators to consider per measurement
# Fourth argument tells us how many smearings to consider per measurement
# Fifth argument tells us the directory to analyze
# Final argument tells us whether to consider the Konishi or SUGRA
if len(sys.argv) < 7:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <# ops> <# smear> <dir> <op>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
ops_per_smear = int(sys.argv[3])
num_smear = int(sys.argv[4])
smear = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
num_ops = ops_per_smear * num_smear
tag = str(sys.argv[5])
op = 'O' + str(sys.argv[6]) + ' '
runtime = -time.time()

# Quick sanity check
if ops_per_smear < 1 or ops_per_smear > 2:
  print "ERROR: Have at most 2 operators per measurement"
  sys.exit(1)
if num_smear > len(smear):
  print "ERROR: Only set up", num_smear, "smearings per measurement"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct list of which configurations have been analyzed
cfgs = []
files = tag + '/mcrg.*'
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
firstfile = tag + '/mcrg.' + str(cfgs[0])
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
tBK = np.zeros_like(tAK)
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
  filename = tag + '/mcrg.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  for line in open(toOpen[0]):
    # Format: OK smear bl Tr(PP) Tr(UU)
    # op = 'OK ' or 'OS ' defined at top
    if line.startswith(op):
      temp = line.split()
      this_smear = int(temp[1])
      N = -1
      for i in range(num_smear):
        if this_smear == int(smear[i]):
          N = i
      if N < 0:
        continue
      bl = int(temp[2])
      dat[ops_per_smear * N][bl] = float(temp[3])
      if ops_per_smear > 1:
        dat[ops_per_smear * N + 1][bl] = float(temp[4])

    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print toOpen[0], "reports two measurements"
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
# Optionally load xi^4 for each operator on each blocking level
# or just set it to unity
# In either case ignore uncertainties, which tend to be small
# Jackknife results for xi^4
xi = np.zeros((num_ops, blmax + 1), dtype = np.float)
xi_file = 'results/xi.dat'
if os.path.isfile(xi_file):
  print "Reading xi^4 from", xi_file
  for line in open(xi_file):
    if line.startswith('# '):
      continue
    else:
      temp = line.split()
      bl = int(temp[0])
      for i in range(num_smear):
        # With the polar field we probably only want xi^2...
        xi[ops_per_smear * i][bl] = np.sqrt(float(temp[1]))
        # With the U.Ubar field we probably want xi^4...
        if ops_per_smear > 1:
          xi[ops_per_smear * i + 1][bl] = float(temp[1])
else:
  print xi_file, "does not exist..."
  sys.exit(1)
  print "Setting xi^4 to unity"
  for i in range(num_smear):
    for j in range(ops_per_smear):
      xi[ops_per_smear * i + j][0] = 1.0
      for bl in range(1, blmax + 1):
        xi[ops_per_smear * i + j][bl] = xi[ops_per_smear * i + j][bl - 1] * 1.0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination,
# and analyze them to determine A(n), B(n) and T(n),
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
print "Analyzing %s with %d blocks of length %d MDTU" \
      % (op.rstrip(), Nblocks, block_size)
print "%d operators in the stability matrix" % num_ops
outfile = open('results/MCRG.dat', 'w')
print >> outfile, "# Analyzing %s with %d blocks of length %d MDTU" \
                  % (op.rstrip(), Nblocks, block_size)
for bl in range(blmax + 1):
  print >> outfile, "xi4 %d" % bl,
  for i in range(num_smear - 1):
    print >> outfile, "%.4g" % xi[i][bl],
  print >> outfile, "%.4g" % xi[num_smear - 1][bl]

print >> outfile, "# %d operators in the stability matrix" % num_ops
for bl in range(blmax):
  ave = np.average(jkDeltaK[bl])
  err = (Nblocks - 1.0) * np.sum((jkDeltaK[bl] - ave)**2) / float(Nblocks)
  print >> outfile, "Delta%s %d %.6g %.4g" \
                    % (str(sys.argv[6]), bl + 1, ave, err)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
