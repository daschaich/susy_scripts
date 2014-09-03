#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Run jackknifed fits of zero-momentum-projected SUGRA correlator
#   C_K(t) = C + A(e^{-M_K t} + e^{-M_K (Nt - t)})
# as well as the finite differences for both Konishi and SUGRA,
#   D(t) = A(e^{-M t} - e^{-M (Nt - t)})
# Need to force positive sinh mass, otherwise strange things can happen

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument is the minimum t to use in the fits (t_max = Nt / 2)
# Fourth argument tells us whether to analyze "corr" or "stout" files
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <tmin> <tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tmin = int(sys.argv[3])
tag = str(sys.argv[4])
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

# Extract Nt from first output file
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nt'):
    Nt = int((line.split())[1])
    break

# Now we can define the fit functions for C_S(t) and D(t)
# The err functions will be minimized via least-squares optimization
#cosh = lambda p, t: p[0] * np.cosh(p[1] * (t - (Nt / 2.0)))
cosh = lambda p, t: p[0] * (np.exp(-p[1] * t) + np.exp(-p[1] * (Nt - t)))
cosherr = lambda p, t, y, err: (cosh(p, t) - y) / err
#sinh = lambda p, t: p[0] * np.sinh(p[1] * (t - (Nt / 2.0)))
sinh = lambda p, t: -p[0] * (np.exp(-p[1] * t) - np.exp(-p[1] * (Nt - t)))
# Need to force positive sinh mass, otherwise strange things can happen
def sinherr(p, t, y, err):
  if p[1] > 0:
    return (sinh(p, t) - y) / err
  else:
    return 1e6
p_in = [0.1, 0.1]               # Order-of-magnitude initial guesses

# We exploit t <--> Nt - t symmetry to only print 0 through Nt / 2
Npts = Nt / 2 + 1 - tmin  # Assume Nt is even
dof = Npts - len(p_in)
if dof <= 0:
  print "ERROR: Not enough data points to fit to (2, 2) rational function"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# There is some funny business with indexing from zero despite tmin>0
# We need to shift things by tmin to deal with it
# K = Konishi, S = SUGRA
# DK and DS are corresponding finite differences
Kdat   = [[] for x in range(Npts)]
KdatSq = [[] for x in range(Npts)]
Sdat   = [[] for x in range(Npts)]
SdatSq = [[] for x in range(Npts)]

DKdat   = [[] for x in range(Npts - 1)]
DKdatSq = [[] for x in range(Npts - 1)]
DSdat   = [[] for x in range(Npts - 1)]
DSdatSq = [[] for x in range(Npts - 1)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK   = [0 for x in range(Npts)]
tKSq = [0 for x in range(Npts)]
tS   = [0 for x in range(Npts)]
tSSq = [0 for x in range(Npts)]
tDK   = [0 for x in range(Npts - 1)]
tDKSq = [0 for x in range(Npts - 1)]
tDS   = [0 for x in range(Npts - 1)]
tDSSq = [0 for x in range(Npts - 1)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    if count == 0:
      print "ERROR: no data to average after file %s:" % toOpen
      sys.exit(1)
    for t in range(Npts):
      Kdat[t].append(tK[t] / float(count))
      KdatSq[t].append(tKSq[t] / float(count))
      Sdat[t].append(tS[t] / float(count))
      SdatSq[t].append(tSSq[t] / float(count))
    for t in range(Npts - 1):
      DKdat[t].append(tDK[t] / float(count))
      DKdatSq[t].append(tDKSq[t] / float(count))
      DSdat[t].append(tDS[t] / float(count))
      DSdatSq[t].append(tDSSq[t] / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK   = [0 for x in range(Npts)]
    tKSq = [0 for x in range(Npts)]
    tS   = [0 for x in range(Npts)]
    tSSq = [0 for x in range(Npts)]
    tDK   = [0 for x in range(Npts - 1)]
    tDKSq = [0 for x in range(Npts - 1)]
    tDS   = [0 for x in range(Npts - 1)]
    tDSSq = [0 for x in range(Npts - 1)]

  # Running averages
  prev_time = float('NaN')    # To make it obvious if this isn't overwritten
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: KONISHI t dat
    if line.startswith('KONISHI '):
      temp = line.split()
      t = int(temp[1])
      dat = float(temp[2])
      if t < tmin:
        continue
      elif t == tmin:
        count += 1    # Only increment once per measurement!
      elif t > tmin:
        tDK[t - tmin - 1] += dat - prev_time
        tDKSq[t - tmin - 1] += (dat - prev_time)**2
      prev_time = dat
      tK[t - tmin] += dat
      tKSq[t - tmin] += dat**2
    # We go through all Konishi data before reaching first SUGRA line
    # Format: SUGRA t dat
    elif line.startswith('SUGRA '):
      temp = line.split()
      t = int(temp[1])
      dat = float(temp[2])
      if t < tmin:
        continue
      elif t > tmin:    # Need to reset prev_time below
        tDS[t - tmin - 1] += dat - prev_time
        tDSSq[t - tmin - 1] += (dat - prev_time)**2
      prev_time = dat
      tS[t - tmin] += dat
      tSSq[t - tmin] += dat**2

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  if count == 0:
    print "ERROR: no data to average after file %s:" % toOpen
    sys.exit(1)
  for t in range(Npts):
    Kdat[t].append(tK[t] / float(count))
    KdatSq[t].append(tKSq[t] / float(count))
    Sdat[t].append(tS[t] / float(count))
    SdatSq[t].append(tSSq[t] / float(count))
  for t in range(1, Npts):
    DKdat[t - 1].append(tDK[t - 1] / float(count))
    DKdatSq[t - 1].append(tDKSq[t - 1] / float(count))
    DSdat[t - 1].append(tDS[t - 1] / float(count))
    DSdatSq[t - 1].append(tDSSq[t - 1] / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
# Require multiple blocks for now
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

times  = np.arange(tmin, Nt / 2 + 1, dtype = np.int)
Ktot   = np.array([sum(Kdat[t]) for t in range(Npts)])
KtotSq = np.array([sum(KdatSq[t]) for t in range(Npts)])
Stot   = np.array([sum(Sdat[t]) for t in range(Npts)])
StotSq = np.array([sum(SdatSq[t]) for t in range(Npts)])

# For now, require that times be evenly spaced
dt = times[1] - times[0]
Dtimes = np.empty(Nt / 2 - tmin, dtype = np.float)
for t in range(len(Dtimes)):
  if (times[t + 1] - times[t]) != dt:
    print "ERROR: requiring evenly spaced times but have",
    print times
    sys.exit(1)
  Dtimes[t] = times[t] + 0.5 * float(dt)
DKtot   = np.array([sum(DKdat[t]) for t in range(Npts - 1)])
DKtotSq = np.array([sum(DKdatSq[t]) for t in range(Npts - 1)])
DStot   = np.array([sum(DSdat[t]) for t in range(Npts - 1)])
DStotSq = np.array([sum(DSdatSq[t]) for t in range(Npts - 1)])

# Monitor how many times the fits do or don't converge
NS = 0
NDK = 0
NDS = 0

# Fit results for all jk samples
Sout = np.empty((len(p_in), Nblocks), dtype = np.float)
DKout = np.empty((len(p_in), Nblocks), dtype = np.float)
DSout = np.empty((len(p_in), Nblocks), dtype = np.float)
for i in range(Nblocks):    # Jackknife samples
  # First cosh fits to SUGRA correlator
  S    = np.empty(Npts, dtype = np.float)
  Serr = np.empty(Npts, dtype = np.float)
  for t in range(Npts):
    S[t] = (Stot[t] - Sdat[t][i]) / (Nblocks - 1.)
    temp = (StotSq[t] - SdatSq[t][i]) / (Nblocks - 1.)
    Serr[t] = np.sqrt((temp - S[t]**2) / (Nblocks - 1.))
  success = -1    # In case fit dies without returning
  temp, cov, infodict, mesg, success = \
              optimize.leastsq(cosherr, p_in[:], args=(times, S, Serr), \
                               full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    for j in range(len(p_in)):
      Sout[j][i] = temp[j]
    NS += 1
  else:                               # Ignore failed fits
    print mesg

  # Now sinh fits to finite differences of Konishi and SUGRA correlators
  DK    = np.empty(Npts - 1, dtype = np.float)
  DKerr = np.empty(Npts - 1, dtype = np.float)
  DS    = np.empty(Npts - 1, dtype = np.float)
  DSerr = np.empty(Npts - 1, dtype = np.float)
  for t in range(Npts - 1):
    DK[t] = (DKtot[t] - DKdat[t][i]) / (Nblocks - 1.)
    temp = (DKtotSq[t] - DKdatSq[t][i]) / (Nblocks - 1.)
    DKerr[t] = np.sqrt((temp - DK[t]**2) / (Nblocks - 1.))
    DS[t] = (DStot[t] - DSdat[t][i]) / (Nblocks - 1.)
    temp = (DStotSq[t] - DSdatSq[t][i]) / (Nblocks - 1.)
    DSerr[t] = np.sqrt((temp - DS[t]**2) / (Nblocks - 1.))

  success = -1
  temp, cov, infodict, mesg, success = \
              optimize.leastsq(sinherr, p_in[:], args=(Dtimes, DK, DKerr), \
                               full_output=True)
  if success in (1, 2, 3, 4):
    for j in range(len(p_in)):
      DKout[j][i] = temp[j]
    NDK += 1
  else:
    print mesg

  success = -1
  temp, cov, infodict, mesg, success = \
              optimize.leastsq(sinherr, p_in[:], args=(Dtimes, DS, DSerr), \
                               full_output=True)
  if success in (1, 2, 3, 4):
    for j in range(len(p_in)):
      DSout[j][i] = temp[j]
    NDS += 1
  else:
    print mesg

  # Optionally monitor chiSq and confidence level of fits
  # This requires using the correctly normalized standard error
#  chiSq = (infodict['fvec']**2).sum()
#  CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
#  print "%.4g %d --> %.4g" % (chiSq, dof, CL)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
# Ignore any failed fits
print "Fitting with %d blocks of length %d MDTU..." % (Nblocks, block_size)
print "%d of %d SUGRA fits succeeded" % (NS, Nblocks)
print "%d of %d Konishi finite difference fits succeeded" % (NDK, Nblocks)
print "%d of %d SUGRA finite difference fits succeeded" % (NDS, Nblocks)

outfile = open('results/fit.corr', 'w')
print >> outfile, "# Fitting with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "%d of %d SUGRA fits succeeded" % (NS, Nblocks)
print >> outfile, "%d of %d Konishi finite difference fits succeeded" \
                  % (NDK, Nblocks)
print >> outfile, "%d of %d SUGRA finite difference fits succeeded" \
                  % (NDS, Nblocks)

tag = ["amplitude", "mass"]
if NS > 0:
  for i in range(len(p_in)):
    ave = np.average(Sout[i])
    var = (NS - 1.) * np.sum((Sout[i] - ave)**2) / float(NS)
    print >> outfile, "SUGRA %s %.6g %.4g" % (tag[i], ave, np.sqrt(var))
if NDK > 0:
  for i in range(len(p_in)):
    ave = np.average(DKout[i])
    var = (NDK - 1.) * np.sum((DKout[i] - ave)**2) / float(NDK)
    print >> outfile, "DKonishi %s %.6g %.4g" % (tag[i], ave, np.sqrt(var))
if NDS > 0:
  for i in range(len(p_in)):
    ave = np.average(DSout[i])
    var = (NDS - 1.) * np.sum((DSout[i] - ave)**2) / float(NDS)
    print >> outfile, "DSUGRA %s %.6g %.4g" % (tag[i], ave, np.sqrt(var))

runtime += time.time()
print >> outfile, "# Runtime: %.2g seconds" % runtime

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
