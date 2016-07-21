#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Blocked fits of Wilson loops to W(r, t) = w * exp(-V(r) * t)
# followed by fits of V(r) to the confining form
#   r * V(r) = A * r + C + sigma * r^2
# Only print out sigma to check whether or not it is zero

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
  outfilename = 'results/sigma.dat'
elif loop == 'D':
  outfilename = 'results/sigma_det.dat'
elif loop == 'POLAR':
  outfilename = 'results/sigma_pol.dat'
else:
  print "ERROR:", loop, "unrecognized"
  print "Must be one of POT, D or POLAR"
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
expfunc = lambda p, x: p[0] * np.exp(-p[1] * x)
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = [0.01, 0.1]   # Order-of-magnitude initial guesses
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

# Cycle through first output file to read MAX_T,
# each scalar displacement and total number of them
temp_r = []
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nt'):
    # Smaller MAX_T seems to increase central values
    # !!! May want to include as systematic effect
    MAX_T = int((line.split())[1]) / 2 - 1
    print "MAX_T =", MAX_T
  # Format: hvy_pot: MAX_T = $MAX_T, MAX_X = $MAX_X --> r < $MAX_r
  elif line.startswith('hvy_pot: MAX_T '):
    temp = line.split()
    max_meas = int((temp[3]).rstrip(','))    # Strip ',' from end
    if max_meas < MAX_T:
      print "ERROR: Only measured Wilson loops to %d but MAX_T=%d" \
            % (max_meas, MAX_T)
      sys.exit(1)
      print

  # Format: $tag_LOOP # r t dat
  elif line.startswith(loop + '_LOOP '):
    temp = line.split()
    if int(temp[3]) > 1:      # Only consider first t=1!
      break
    temp_r.append(float(temp[2]))
Npts = len(temp_r)
x_r = np.array(temp_r, dtype = np.float)
#print x_r
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
        Wdat[j][t].append(tW[j][t] / float(count))
        WdatSq[j][t].append(tWSq[j][t] / float(count))
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
  check = -1
  for line in open(toOpen[0]):
    # Format: $tag_LOOP # r t dat
    if line.startswith(loop + '_LOOP '):
      temp = line.split()
      j = int(temp[1])
      t = int(temp[3]) - 1      # Shift to index from zero
      dat = float(temp[4])
      if j == 0 and t == 0:
        count += 1              # Only increment once per measurement!
      if t < MAX_T:
        tW[j][t] += dat
        tWSq[j][t] += dat**2
    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print toOpen[0], "reports two measurements"
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for j in range(Npts):
    for t in range(MAX_T):
      Wdat[j][t].append(tW[j][t] / float(count))
      WdatSq[j][t].append(tWSq[j][t] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Wdat[0][0])

# Sort r, Wdat and WdatSq to reproduce old results
# There ought to be a less awkward way to do this
order = x_r.argsort()
x_r = x_r[order]

tmp = np.zeros(Npts, dtype = np.float)
tmpSq = np.zeros_like(tmp)
for t in range(MAX_T):
  for i in range(Nblocks):
    for j in range(Npts):
      tmp[j] = Wdat[j][t][i]
      tmpSq[j] = WdatSq[j][t][i]
    for j in range(Npts):
      Wdat[j][t][i] = tmp[order[j]]
      WdatSq[j][t][i] = tmpSq[order[j]]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination,
# and fit each Wilson loop to W(r, t) = w * exp(-V(r) * t)
# and then fit r * V(r) = A * r + C + sigma * r^2 to find the string tension
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

# All jackknife results
Nt = MAX_T - 3    # How many t_min we can consider
jkC = np.zeros((Nt, Nblocks), dtype = np.float)
jkA = np.zeros_like(jkC)
jksigma = np.zeros_like(jkC)

# All fits
for t_min in range(1, MAX_T - 2):         # Doesn't include MAX_T - 2
  index = t_min - 1
  x_t = np.arange(t_min, MAX_T + 1)       # At least four points
#  print x_t

  for i in range(Nblocks):  # Jackknife samples
    Vlist = [[] for x in range(Npts)]
    rV = np.empty(Npts, dtype = np.float)
    weight = np.empty_like(rV)

    # Inner jackknife:
    # Fit W(r, t) = w * exp(-V(r) * t) for t_min <= t <= MAX_T
    # Recall that t is indexed from zero instead of one
    Ninner = Nblocks - 1
    WIn = np.empty_like(Wtot)
    WInSq = np.empty_like(Wtot)
    for j in range(Npts):
      for t in range(MAX_T):
        WIn[j][t] = Wtot[j][t] - Wdat[j][t][i]
        WInSq[j][t] = WtotSq[j][t] - WdatSq[j][t][i]
    for ii in range(Nblocks):
      if ii == i:
        continue
      for j in range(Npts):
        W = np.empty(len(x_t), dtype = np.float)
        WErr = np.empty_like(W)
        for t in range(t_min, MAX_T + 1):
          W[t - t_min] = (WIn[j][t - 1] - Wdat[j][t - 1][ii]) / (Ninner - 1.0)
          temp = (WInSq[j][t - 1] - WdatSq[j][t - 1][ii]) / (Ninner - 1.0)
          WErr[t - t_min] = np.sqrt(temp - W[t - t_min]**2)

        # V(r) is second of two parameters returned as temp
        temp, success = optimize.leastsq(errfunc, p_in[:], args=(x_t, W, WErr))
        Vlist[j].append(float(temp[1]))

    # Now we have an estimate of V and weight
    if not len(Vlist[0]) == Ninner:
      print "ERROR: Wrong number of samples in inner jackknife", Ninner,
      print "vs.", len(Vlist[0])
    for j in range(Npts):
      temp = np.array(Vlist[j])
      ave = np.average(temp)
      rV[j] = x_r[j] * ave
      var = (Ninner - 1.0) * np.sum((temp - ave)**2) / float(Ninner)
      weight[j] = 1.0 / (x_r[j] * np.sqrt(var))    # Squared in fit
#     print x_r[j], rV[j][index][i], weight[j]

    # Fit r * V(r) = sigma * r^2 + A * r - C for all r
    # [sigma, A, -C] are the three parameters returned as temp
    temp = np.polyfit(x_r, rV, 2, full=False, w=weight)
    jksigma[index][i] = temp[0]
    jkA[index][i] = temp[1]
    jkC[index][i] = -1.0 * temp[2]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Analyzing with %d blocks of length %d MDTU" % (Nblocks, block_size)

# Print string tension sigma for all t_min
# Include Coulomb coefficient to see how it changes
# Include asymptotic value (without uncertainty) to check fits offline
outfile = open(outfilename, 'w')
print >> outfile, "# Analyzing with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)

for t_min in range(Nt):
  ave = np.average(jksigma[t_min])
  var = (Nblocks - 1.0) * np.sum((jksigma[t_min] - ave)**2) / float(Nblocks)
  print >> outfile, "%d %.6g %.4g" % (t_min + 1, ave, np.sqrt(var)),

  ave = np.average(jkC[t_min])
  var = (Nblocks - 1.0) * np.sum((jkC[t_min] - ave)**2) / float(Nblocks)
  print >> outfile, "# %.6g %.4g" % (ave, np.sqrt(var)),
  print >> outfile, "# %.6g " % np.average(jkA[t_min])
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
