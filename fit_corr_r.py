#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Run jackknifed fits of radial correlators
#   C(r) = A' * r^{Delta}

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third and fourth arguments are the fit range (r_min, r_max)
# Fifth argument tells us which files to analyze
#   (for example, "corr", "konishi0", "subtracted0")
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <rmin> <rmax> <tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
rMin = float(sys.argv[3])
rMax = float(sys.argv[4])
tag = str(sys.argv[5])
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
  if cfg not in cfgs and cfg >= cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files", files, "found"
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]

# Cycle through first output file to count how many scalar distances
# we will consider, and record their order
r = []          # List of r
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  # Format: CORR_K label r tag1 tag2 dat_vev dat_vol
  if line.startswith('CORR_K '):
    temp = line.split()
    tag1 = int(temp[3])
    tag2 = int(temp[4])
    if tag1 == 0 and tag2 == 0:   # Only count each r once
      tr = float(temp[2])
      if tr < rMin or tr > rMax:
        continue
      r.append(tr)
  elif line.startswith('CORR_S '):
    break             # Don't go through whole file yet

Npts = len(r)

# Now we can define the fit function
fitfunc = lambda p, x: p[0] / np.power(x, p[1])
errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
p_in = [1.0, 1.0]    # Order-of-magnitude initial guesses
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
# For now only consider log-polar operator with ensemble subtraction
Kdat   = [[] for x in range(Npts)]
KdatSq = [[] for x in range(Npts)]
Sdat   = [[] for x in range(Npts)]
SdatSq = [[] for x in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK   = [0.0 for x in range(Npts)]
tKSq = [0.0 for x in range(Npts)]
tS   = [0.0 for x in range(Npts)]
tSSq = [0.0 for x in range(Npts)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for i in range(Npts):
      Kdat[i].append(tK[i] / float(count))
      KdatSq[i].append(tKSq[i] / float(count))
      Sdat[i].append(tS[i] / float(count))
      SdatSq[i].append(tSSq[i] / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK   = [0.0 for x in range(Npts)]
    tKSq = [0.0 for x in range(Npts)]
    tS   = [0.0 for x in range(Npts)]
    tSSq = [0.0 for x in range(Npts)]

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  rCount_K = 0
  rCount_S = 0
  for line in open(toOpen[0]):
    # Format: CORR_? label r tag1 tag2 dat_vev dat_vol
    if line.startswith('CORR'):
      if line.startswith('CORR_K 0 0 0 0 '):
        count += 1            # Only increment once per measurement!

      # For now only consider log-polar operator with ensemble subtraction
      temp = line.split()
      tr = float(temp[2])
      tag1 = int(temp[3])
      tag2 = int(temp[4])
      if tr < rMin or tr > rMax or not tag1 == 0 or not tag2 == 0:
        continue
      dat = float(temp[5])

      if line.startswith('CORR_K '):
        tK[rCount_K] += dat
        tKSq[rCount_K] += dat**2
        rCount_K += 1
      elif line.startswith('CORR_S '):
        tS[rCount_S] += dat
        tSSq[rCount_S] += dat**2
        rCount_S += 1

    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print infile, "reports two measurements"
      check = 1
      if not rCount_K == Npts or not rCount_S == Npts:
        print infile, "gives wrong number of points"
        sys.exit(1)
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)

# Check special cases
# 1) Last block is full (ssume last few measurements are equally spaced)
# 2) Only a single block
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]) or len(Kdat[0]) < 1:
  for i in range(Npts):
    Kdat[i].append(tK[i] / float(count))
    KdatSq[i].append(tKSq[i] / float(count))
    Sdat[i].append(tS[i] / float(count))
    SdatSq[i].append(tSSq[i] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# If only one block, just do a single fit with propagated uncertainties
# from the covariance matrix
# In this case 'count' remains the total number of measurements
# optimize.leastsq does not scale by chiSq_dof, which we desire
if Nblocks == 1:
  outfile = open('results/fit.corr', 'w')
  print "Single fit with range %.4g <= r <= %.4g" % (rMin, rMax)
  print >> outfile, "Single fit with range %.4g <= r <= %.4g" % (rMin, rMax)

  pts = np.array(r)
  K = np.array([Kdat[x][0] for x in range(Npts)])
  KErr = np.empty_like(K)
  for x in range(Npts):
    err = KdatSq[x][0] - Kdat[x][0]**2
    KErr[x] = np.sqrt(err / (float(count) - 1.0))

  temp, cov, infodict, mesg, success = \
              optimize.leastsq(errfunc, p_in[:], args=(pts, K, KErr),
                               full_output=True)

  Delta = 0.5 * temp[1]
  err = 0.5 * np.sqrt(cov[1][1])
  print >> outfile, "Konishi %.6g %.4g" % (Delta, err)

  S = np.array([Sdat[x][0] for x in range(Npts)])
  SErr = np.empty_like(S)
  for x in range(Npts):
    err = SdatSq[x][0] - Sdat[x][0]**2
    SErr[x] = np.sqrt(err / (float(count) - 1.0))

  temp, cov, infodict, mesg, success = \
              optimize.leastsq(errfunc, p_in[:], args=(pts, S, SErr),
                               full_output=True)

  Delta = 0.5 * temp[1]
  err = 0.5 * np.sqrt(cov[1][1])
  print >> outfile, "SUGRA %.6g %.4g" % (Delta, err)

  runtime += time.time()
  print "# Runtime: %.2g seconds" % runtime
  print >> outfile, "# Runtime: %.2g seconds" % runtime
  outfile.close()
  sys.exit(0)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
columns = [('r', float), ('dat', float), ('err', float)]
Ktot    = np.array([sum(Kdat[x]) for x in range(Npts)])
KtotSq  = np.array([sum(KdatSq[x]) for x in range(Npts)])
Stot    = np.array([sum(Sdat[x]) for x in range(Npts)])
StotSq  = np.array([sum(SdatSq[x]) for x in range(Npts)])

# Fit results for all jk samples
Kout = []
Sout = []
for i in range(Nblocks):    # Jackknife samples
  # Sort data in case that may help fitter
  K = np.zeros(Npts, dtype = columns)
  S = np.zeros_like(K)
  for x in range(Npts):
    K[x][0] = r[x]
    K[x][1] = (Ktot[x] - Kdat[x][i]) / (Nblocks - 1.0)
    temp = (KtotSq[x] - KdatSq[x][i]) / (Nblocks - 1.0)
    K[x][2] = np.sqrt((temp - K[x][1]**2) / (Nblocks - 1.0))

    S[x][0] = r[x]
    S[x][1] = (Stot[x] - Sdat[x][i]) / (Nblocks - 1.0)
    temp = (StotSq[x] - SdatSq[x][i]) / (Nblocks - 1.0)
    S[x][2] = np.sqrt((temp - S[x][1]**2) / (Nblocks - 1.0))

  # Sort -- note that r itself is not sorted
  K = np.sort(K, order='r')
  S = np.sort(S, order='r')

  # Copy Konishi data into appropriate structures and fit
  pts = np.zeros(Npts, dtype = float)
  dat = np.zeros(Npts, dtype = float)
  err = np.zeros(Npts, dtype = float)
  for x in range(Npts):
    pts[x] = K[x][0]
    dat[x] = K[x][1]
    err[x] = K[x][2]
  temp, cov, infodict, mesg, success = \
              optimize.leastsq(errfunc, p_in[:], args=(pts, dat, err), \
                               full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    Kout.append(0.5 * temp[1])
  else:                               # Ignore failed fits
    print mesg

  # Copy SUGRA data into appropriate structures and fit
  for x in range(Npts):
    dat[x] = S[x][1]
    err[x] = S[x][2]
  temp, cov, infodict, mesg, success = \
              optimize.leastsq(errfunc, p_in[:], args=(pts, dat, err), \
                               full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    Sout.append(0.5 * temp[1])
  else:                               # Ignore failed fits
    print mesg

  # The uncorrelated chiSq isn't meaningful
  # but might be worth monitoring to make sure it is <<1
#  chiSq = (infodict['fvec']**2).sum()
#  print "chiSq %.4g" % chiSq
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
# Ignore any failed fits
Kpow = np.array(Kout)
Spow = np.array(Sout)
NK = len(Kpow)
NS = len(Spow)

print "Fitting with %d blocks of length %d MDTU..." % (Nblocks, block_size)
print "Fit range %.4g <= r <= %.4g contains %d points" % (rMin, rMax, Npts)
print "%d of %d Konishi fits succeeded" % (NK, Nblocks)
print "%d of %d SUGRA fits succeeded" % (NS, Nblocks)

outfile = open('results/fit.corr', 'w')
print >> outfile, "# Fitting with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "Fit range %.4g <= r <= %.4g contains %d points" \
                  % (rMin, rMax, Npts)
print >> outfile, "%d of %d Konishi fits succeeded" % (NK, Nblocks)
print >> outfile, "%d of %d SUGRA fits succeeded" % (NS, Nblocks)

if NK > 0:
  ave = np.average(Kpow)
  var = (NK - 1.0) * np.sum((Kpow - ave)**2) / float(NK)
  print >> outfile, "Konishi %.6g %.4g" % (ave, np.sqrt(var))
if NS > 0:
  ave = np.average(Spow)
  var = (NS - 1.0) * np.sum((Spow - ave)**2) / float(NS)
  print >> outfile, "SUGRA %.6g %.4g" % (ave, np.sqrt(var))

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime
print >> outfile, "# Runtime: %.2g seconds" % runtime

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
