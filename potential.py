#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy.optimize import least_squares
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

# errfunc will be minimized via least-squares optimization
# p_in are order-of-magnitude initial guesses
expfunc = lambda p, x: p[0] * np.exp(-p[1] * x)
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = np.array([0.01, 0.1])

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  J[:, 0] = np.exp(-p[1] * x)
  J[:, 1] = -p[0] * x * np.exp(-p[1] * x)
  for i in range(p.size):
    J[:, i] /= err
  return J

# Another least squares setup for the Coulomb potential fit
#   r * V(r) = A * r - C
# This is overkill for a linear fit
# but as of November 2018 np.polyfit has a covariance issue
# github.com/numpy/numpy/issues/11196
linfunc = lambda V, x: V[0] * x - V[1]
err_lin = lambda V, x, y, err: (linfunc(V, x) - y) / err
V_in = np.array([0.01, 0.01])

# Define corresponding Jacobian matrix
def jac_lin(V, x, y, err):
  J = np.empty((x.size, V.size), dtype = np.float)
  J[:, 0] = x
  J[:, 1] = -1.0
  for i in range(V.size):
    J[:, i] /= err
  return J
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
# and then fit r * V(r) = A * r - C to find the Coulomb coefficients
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
Nt = MAX_T - 2    # How many t_min we can consider
jkV = np.zeros((Npts, Nt, Nblocks), dtype = np.float)
jkC = np.zeros((Nt, Nblocks), dtype = np.float)
jkA = np.zeros_like(jkC)

# !!! Temporary hack to provide jackknife ratios
# Will write each jackknife estimate to this file for post-processing
tempfilename = 'results/C' + loop + '.jk'
tempfile = open(tempfilename, 'w')

# All fits
for t_min in range(1, MAX_T - 1):         # Doesn't include MAX_T - 1
  index = t_min - 1
  x_t = np.arange(t_min, MAX_T + 1)       # At least three points
#  print x_t

  for i in range(Nblocks):  # Jackknife samples
    jkVlist = [[] for x in range(Npts)]
    rV = np.empty(Npts, dtype = np.float)
    rVerr = np.empty_like(rV)

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

        # V(r) is second of two parameters returned as all_out.x
        # Simply ignore covariance matrix at this stage
        # method='lm' is Levenberg--Marquardt (can't handle bounds)
        all_out = least_squares(errfunc, p_in, jac=jac, max_nfev=10000,
                                method='lm', args=(x_t, W, WErr))
        temp = all_out.x
        if all_out.success < 0 or all_out.success > 4:
          print("ERROR: Fit failed with the following error message")
          print(errmsg)
          sys.exit(1)

        jkVlist[j].append(float(temp[1]))

    # Now we have an estimate of jkV, rV and rVerr
    if not len(jkVlist[0]) == Ninner:
      print "ERROR: Wrong number of samples in inner jackknife", Ninner,
      print "vs.", len(jkVlist[0])
    for j in range(Npts):
      temp = np.array(jkVlist[j])
      jkV[j][index][i] = np.average(temp)
      rV[j] = x_r[j] * jkV[j][index][i]
      var = (Ninner - 1.0) * np.sum((temp - jkV[j][index][i])**2) / float(Ninner)
      rVerr[j] = x_r[j] * np.sqrt(var)
#     print x_r[j], jkV[j][index][i], rVerr[j]

    # Fit r * V(r) = A * r - C for all r
    # [A, C] are the two parameters returned as all_out.x
    all_out = least_squares(err_lin, V_in, jac=jac_lin, max_nfev=10000,
                            method='lm', args=(x_r, rV, rVerr))
    temp = all_out.x
    if all_out.success < 0 or all_out.success > 4:
      print("ERROR: Fit failed with the following error message")
      print(errmsg)
      sys.exit(1)

    jkA[index][i] = temp[0]
    jkC[index][i] = temp[1]

    # !!! Temporary hack to provide jackknife ratios
    if MAX_T == 11 and t_min == 7:
      print >> tempfile, "%.6g" % jkC[index][i]
    elif MAX_T == 15 and t_min == 12:
      print >> tempfile, "%.6g" % jkC[index][i]
tempfile.close()
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
  print >> Vfile, "%.4g" % x_r[j],
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
  print >> outfile, "%.4g err            " % x_r[j],
print >> outfile, "%.4g err" % x_r[Npts - 1]

for t in range(MAX_T):
  print >> outfile, "%d" % (t + 1),     # Shift since indexed from zero
  for j in range(Npts):
    dat = np.array(Wdat[j][t])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
    print >> outfile, "%.6g %.4g" % (ave, err),
  print >> outfile, ""
outfile.close()

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
