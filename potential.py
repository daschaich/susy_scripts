#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Blocked fits of Wilson loops to W(r, t) = w * exp(-V(r) * t)
# followed by fits of V(r) to the Coulomb form r * V(r) = A * r + C
# Print out Coulomb coefficients C and every V(r) as functions of
# t_min in the fits, as well as average W(r, t) themselves

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us whether to analyze "corr" or "stout" files
# Fourth argument tells us which type of loops to consider: POT, D or POLAR
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <tag> <loop>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
loop = str(sys.argv[4])
runtime = -time.time()

# Check loop tag
if loop == 'POT':
  Cfilename = 'results/C.dat'
  Vfilename = 'results/V.dat'
  outfilename = 'results/pot_loops.dat'
elif loop == 'D':
  Cfilename = 'results/Cdet.dat'
  Vfilename = 'results/Vdet.dat'
  outfilename = 'results/det_loops.dat'
elif loop == 'POLAR':
  Cfilename = 'results/Cpol.dat'
  Vfilename = 'results/Vpol.dat'
  outfilename = 'results/polar_loops.dat'
else:
  print "ERROR:", loop, "unrecognized"
  print "Must be one of POT, D or POLAR"
  sys.exit(1)

# Convenience constants
TOL = 1.0e-6
invSq2  = 1.0 / np.sqrt(2)
invSq6  = 1.0 / np.sqrt(6)
invSq12 = 1.0 / np.sqrt(12)

# errfunc will be minimized via least-squares optimization
expfunc = lambda p, x: p[0] * np.exp(-p[1] * x)
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = [0.01, 0.1]   # Order-of-magnitude initial guesses
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Map (x, y, z) to r on L^3 reduction of A4* lattice
# Check all possible periodic shifts to find true r
def A4map(x_in, y_in, z_in, L):
  r = 100.0 * L  # To be overwritten
  for x in [x_in + L, x_in, x_in - L]:
    for y in [y_in + L, y_in, y_in - L]:
      for z in [z_in + L, z_in, z_in - L]:
        test = np.sqrt((x**2 + y**2 + z**2) * 0.75 \
                       - (x * (y + z) + y * z) * 0.5)

        # Sanity check -- can be commented out for more speed
#        x_a4 = (x - y) * invSq2
#        y_a4 = (x + y - 2.0 * z) * invSq6
#        z_a4 = (x + y + z) * invSq12
#        check = np.sqrt(x_a4**2 + y_a4**2 + z_a4**2)
#        if np.fabs(test - check) > TOL:
#          print "ERROR: %.4g isn't %.4g for (%d, %d, %d)" \
#                % (check, test, x, y, z)
#          sys.exit(1)

        # True r is the smallest
        if test < r:
          r = test
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
# and determine their multiplicity, also recording MAX_T
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
  # Format: hvy_pot: MAX_T = #, MAX_X = #
  elif line.startswith('hvy_pot: MAX_T '):
    temp = line.split()
    MAX_T = int((temp[3]).rstrip(','))    # Strip ',' from end
    MAX_X = int(temp[6])
    num_y = 2 * MAX_X + 1
    lookup = np.empty((MAX_X + 1, num_y, num_y), dtype = np.float)

    # Decide where to cut r
    MAX_r = 100.0 * MAX_X   # To be overwritten
    for y in range(MAX_X + 1):
      for z in range(MAX_X + 1):
        test = A4map(MAX_X + 1, y, z, L)
        if test < MAX_r:
          MAX_r = test
#    print MAX_r

  # Format: POT_LOOP x y z t dat      (similarly for D_LOOP and POLAR_LOOP)
  elif line.startswith(loop + '_LOOP '):
    temp = line.split()
    if int(temp[4]) > 1:
      break                 # Only consider first t=1!
    x = int(temp[1])
    y = int(temp[2])
    z = int(temp[3])
    this_r = A4map(x, y, z, L)

    # Save this_r in lookup table so we can avoid calling A4map in the future
    lookup[x][y + MAX_X][z + MAX_X] = this_r

    # Accumulate multiplicities if 0.5 < this_r < MAX_r
    # Try to avoid roundoff issues in latter comparison
    if this_r < 0.5 or this_r - MAX_r > -TOL:
      continue
    done = -1
    for j in range(len(r)):
      if np.fabs(r[j][0] - this_r) < TOL:
        r[j][1] += 1
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
# Construct arrays of blocked measurements for specified loop
Wdat = [[[] for i in range(MAX_T)] for i in range(Npts)]
WdatSq = [[[] for i in range(MAX_T)] for i in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tW = np.zeros((Npts, MAX_T), dtype = np.float)
tWSq = np.zeros_like(tW)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for j in range(Npts):
      for t in range(MAX_T):
        Wdat[j][t].append(tW[j][t] / float(count * r[j][1]))
        WdatSq[j][t].append(tWSq[j][t] / float(count * r[j][1]))
        tW[j][t] = 0.0
        tWSq[j][t] = 0.0
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
  for line in open(toOpen[0]):
    # Format: *_LOOP x y z t dat
    if line.startswith(loop + '_LOOP '):
      temp = line.split()
      x = int(temp[1])
      y = int(temp[2])
      z = int(temp[3])
      t = int(temp[4]) - 1    # Shift to index from zero
      dat = float(temp[5])
      if x == 0 and y == 0 and z == 0 and t == 0:
        count += 1            # Only increment once per measurement!
      this_r = lookup[x][y + MAX_X][z + MAX_X]

      # Accumulate data if this_r < MAX_r
      # Try to avoid roundoff shenanigans
      if this_r < 0.5 or this_r - MAX_r > -TOL:
        continue
      done = -1
      for j in range(len(r)):
        if np.fabs(r[j][0] - this_r) < TOL:
          done = 1
          tW[j][t] += dat
          tWSq[j][t] += dat**2
          break
      if done < 0:
        print "ERROR: displacement", this_r, "not found"
        sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for j in range(Npts):
    for t in range(MAX_T):
      Wdat[j][t].append(tW[j][t] / float(count * r[j][1]))
      WdatSq[j][t].append(tWSq[j][t] / float(count * r[j][1]))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Wdat[0][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination,
# and fit each Wilson loop to W(r, t) = w * exp(-V(r) * t)
# and then fit r * V(r) = A * r + C to find the Coulomb coefficients
# Require multiple blocks instead of worrying about error propagation
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

Wtot = np.empty_like(tW)
WtotSq = np.empty_like(tW)
for j in range(Npts):
  for t in range(MAX_T):
    Wtot[j][t] = sum(Wdat[j][t])
    WtotSq[j][t] = sum(WdatSq[j][t])

x_r = np.empty(Npts, dtype = np.float)
for j in range(Npts):
  x_r[j] = r[j][0]
#print x_r

# All jackknife results
Nt = MAX_T - 2    # How many t_min we can consider
jkV = np.zeros((Npts, Nt, Nblocks), dtype = np.float)
jkC = np.zeros((Nt, Nblocks), dtype = np.float)
jkA = np.zeros_like(jkC)

# !!! Temporary hack to provide jackknife ratios
# Will write each jackknife estimate to this file for post-processing
#tempfilename = 'data/C' + loop + '.jk'
#tempfile = open(tempfilename, 'w')

# All fits
for t_min in range(1, MAX_T - 1):         # Doesn't include MAX_T - 1
  index = t_min - 1
  x_t = np.arange(t_min, MAX_T + 1)       # At least three points
#  print x_t

  for i in range(Nblocks):  # Jackknife samples
    rV = np.empty(Npts, dtype = np.float)
    weight = np.empty_like(rV)

    # Fit W(r, t) = w * exp(-V(r) * t) for t_min <= t <= MAX_T
    # Recall that t is indexed from zero instead of one
    for j in range(Npts):
      W = np.empty(len(x_t), dtype = np.float)
      WErr = np.empty_like(W)
      for t in range(t_min, MAX_T + 1):
        W[t - t_min] = (Wtot[j][t - 1] - Wdat[j][t - 1][i]) / (Nblocks - 1.0)
        temp = (WtotSq[j][t - 1] - WdatSq[j][t - 1][i]) / (Nblocks - 1.0)
        WErr[t - t_min] = np.sqrt(temp - W[t - t_min]**2)

      # V(r) is second of two parameters returned as temp[0]
      temp = optimize.leastsq(errfunc, p_in[:], args=(x_t, W, WErr),
                              full_output = 1)
      jkV[j][index][i] = temp[0][1]
      rV[j] = r[j][0] * jkV[j][index][i]
      weight[j] = 1.0 / (r[j][0] * np.sqrt(temp[1][1][1]))   # Squared in fit...
#      print x_r[j], jkV[j][index][i], weight[j]

    # Fit r * V(r) = A * r - C for all r
    # [A, -C] are the two parameters returned as temp
    temp = np.polyfit(x_r, rV, 1, full=False, w=weight, cov=False)
    jkA[index][i] = temp[0]
    jkC[index][i] = -1.0 * temp[1]

    # !!! Temporary hack to provide jackknife ratios
#    if t_min == 6:
#      print >> tempfile, "%.6g" % jkC[index][i]
#tempfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)

# Start with Coulomb coefficients for all t_min
# Also record asymptotic value to check fits (don't care about its error)
Cfile = open(Cfilename, 'w')
print >> Cfile, "# Analyzing with %d blocks of length %d MDTU" \
                % (Nblocks, block_size)
for t_min in range(Nt):
  ave = np.average(jkC[t_min])
  var = (Nblocks - 1.0) * np.sum((jkC[t_min] - ave)**2) / float(Nblocks)
  print >> Cfile, "%d %.6g %.4g" % (t_min + 1, ave, np.sqrt(var)),

  ave = np.average(jkA[t_min])
  print >> Cfile, "# %.4g" % ave
Cfile.close()

# Now record results for potential V(r) vs. r
# Put each t_min in a separate column
Vfile = open(Vfilename, 'w')
print >> Vfile, "# Analyzing with %d blocks of length %d MDTU" \
                % (Nblocks, block_size)
print >> Vfile, "# r",
for t_min in range(1, MAX_T - 2):
  print >> Vfile, "%d err            " % t_min,
print >> Vfile, "%d err" % (MAX_T - 2)

for j in range(Npts):
  print >> Vfile, "%.4g" % r[j][0],
  for t in range(Nt):
    ave = np.average(jkV[j][t])
    var = (Nblocks - 1.0) * np.sum((jkV[j][t] - ave)**2) / float(Nblocks)
    print >> Vfile, "%.6g %.4g" % (ave, np.sqrt(var)),
  print >> Vfile, ""
Vfile.close()

# Finally save average values of all the loops W(r, t) vs. t
# Put each r in a separate column
outfile = open(outfilename, 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "# t",
for j in range(Npts - 1):
  print >> outfile, "%.4g err            " % r[j][0],
print >> outfile, "%.4g err" % r[Npts - 1][0]

for t in range(MAX_T):
  print >> outfile, "%d" % (t + 1),     # Shift since indexed from zero
  for j in range(Npts):
    dat = np.array(Wdat[j][t])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
    print >> outfile, "%.6g %.4g" % (ave, err),
  print >> outfile, ""
outfile.close()

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
runtime += time.time()
print "Runtime: %.2g seconds" % runtime
# ------------------------------------------------------------------
