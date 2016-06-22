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
#   C(r) = A * r^{Delta}

# For now only consider log-polar scalar field

# Parse arguments: first two are the fit range (r_min, r_max)
# The final argument tells us which file to analyze
# This file already handles the thermalization cut and blocking
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<rmin> <rmax> <file>"
  sys.exit(1)
rMin = float(sys.argv[1])
rMax = float(sys.argv[2])
toOpen = str(sys.argv[3])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Define the fit function
fitfunc = lambda p, x: p[0] / np.power(x, p[1])
errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
p_in = [1.0, 1.0]    # Order-of-magnitude initial guesses

# Take a first pass through the file to read the number of blocks,
# count how many scalar distances we will consider, and record their order
r = []          # List of r
for line in open(toOpen):
  if line.startswith('Nblock '):
    Nblocks = int((line.split())[1])

  # !!! Assume measurements always separated by 10 MDTU
  elif line.startswith('Nmeas '):
    block_size = 10 * int((line.split())[1])

  # Format: BLOCK_K block# label r tag1 tag2 dat
  elif line.startswith('BLOCK_K 0 '):
    temp = line.split()
    tag1 = int(temp[4])
    tag2 = int(temp[5])
    if tag1 == 0 and tag2 == 0:   # Only count each r once
      tr = float(temp[3])
      if tr < rMin or tr > rMax:
        continue
      r.append(tr)
  elif line.startswith('BLOCK_S 0 '):
    break             # Don't go through whole file yet

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read in blocked measurements of each correlator
# K = Konishi, S = SUGRA (averaged over all independent components)
# For now only consider log-polar scalar field
Kdat   = np.empty((Npts, Nblocks), dtype = np.float)
KdatSq = np.empty_like(Kdat)
Sdat   = np.empty_like(Kdat)
SdatSq = np.empty_like(Kdat)

check = -1
oldBlock = 0
rCount_K = 0
rCount_S = 0
for line in open(toOpen):
  # Format: BLOCK_K block# label r tag1 tag2 dat
  if line.startswith('BLOCK'):
    temp = line.split()
    block = int(temp[1])
    tr = float(temp[3])
    tag1 = int(temp[4])
    tag2 = int(temp[5])
    if tr < rMin or tr > rMax or not tag1 == 0 or not tag2 == 0:
      continue
    dat = float(temp[6])

    # Check and reset counters if we've moved on to the next block
    if not block == oldBlock:
      if not rCount_K == Npts or not rCount_S == Npts:
        print "Wrong number of points in block", oldBlock
      oldBlock = block
      rCount_K = 0
      rCount_S = 0

    if line.startswith('BLOCK_K '):
      Kdat[rCount_K][block] = dat
      KdatSq[rCount_K][block] = dat**2
      rCount_K += 1
    elif line.startswith('BLOCK_S '):
      Sdat[rCount_S][block] = dat
      SdatSq[rCount_S][block] = dat**2
      rCount_S += 1

  elif line.startswith('RUNNING COMPLETED'):
    if check == 1:    # Check that we have one measurement per file
      print infile, "reports two measurements"
    check = 1
if check == -1:
  print toOpen, "did not complete"
  sys.exit(1)
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
outfile.close()
# ------------------------------------------------------------------
