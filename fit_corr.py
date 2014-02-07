#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Run jackknifed cosh fits of Konishi and SUGRA correlators

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct list of which configurations have been analyzed
cfgs = []
for filename in glob.glob('Out/corr.*'):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs and cfg >= cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files Out/corr.* found"
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]

# Extract lattice volume from path
# For now only Nt is used; we assume it is a two-digit number!!!
path = os.getcwd()
if 'Gauge' in path:   # Tom's runs
  temp = path.split('Gauge')
  L = int(temp[1][:1])      # First digit after 'Gauge'
  Nt = int(temp[1][1:3])    # Second and third digits after 'nt'
else:                 # My runs
  temp = path.split('nt')
  L = int((temp[0].split('_'))[-1])     # Everything between '_' and 'nt'
  Nt = int(temp[1][:2])                 # First two digits after 'nt'

# Now we can define the constant+cosh around Nt / 2
# (We have some funny business with dropping t=0 but indexing from zero...)
# errfunc will be minimized via least-squares optimization
Kcosh = lambda p, t: p[0] + p[1] * np.cosh(p[2] * (t + 1 - Nt / 2))
Kerrfunc = lambda p, t, y, err: (Kcosh(p, t) - y) / err
K_in = [20., 0.0001, 1.]          # Order-of-magnitude initial guesses
Scosh = lambda p, t: p[0] * np.cosh(p[1] * (t + 1 - Nt / 2))
Serrfunc = lambda p, t, y, err: (Scosh(p, t) - y) / err
S_in = [0.00001, 0.5]             # Order-of-magnitude initial guesses

# We exploit t <--> Nt - t symmetry to only print 0 through Nt / 2
Npts = Nt / 2 + 1  # Assume Nt is even
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# Ignore t=0, which is always well off from any curve
# (However, array will start at zero, so we shift things when reading in...)
# For now just average over all 25 SUGRA components
# K = Konishi, S = SUGRA
Kdat = [[] for x in range(Npts - 1)]
KdatSq = [[] for x in range(Npts - 1)]
Sdat = [[] for x in range(Npts - 1)]
SdatSq = [[] for x in range(Npts - 1)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [0 for x in range(Npts - 1)]
tKSq = [0 for x in range(Npts - 1)]
tS = [0 for x in range(Npts - 1)]
tSSq = [0 for x in range(Npts - 1)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for t in range(Npts - 1):
      Kdat[t].append(tK[t] / float(count))
      KdatSq[t].append(tKSq[t] / float(count))
      Sdat[t].append(tS[t] / float(25. * count))      # Average over all
      SdatSq[t].append(tSSq[t] / float(25. * count))  # Average over all
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK = [0 for x in range(Npts - 1)]
    tKSq = [0 for x in range(Npts - 1)]
    tS = [0 for x in range(Npts - 1)]
    tSSq = [0 for x in range(Npts - 1)]

  # Running averages
  filename = 'Out/corr.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: KONISHI t dat
    if line.startswith('KONISHI '):
      temp = line.split()
      t = int(temp[1])
      if t == 0:
        count += 1    # Only increment once per measurement!
        continue      # Ignore t=0
      dat = float(temp[2])
      tK[t - 1] += dat
      tKSq[t - 1] += dat**2
    # Format: SUGRA a b t dat
    elif line.startswith('SUGRA '):
      temp = line.split()
      t = int(temp[3])
      if t == 0:      # Already incremented count above
        continue      # Ignore t=0
      dat = float(temp[4])
      tS[t - 1] += dat
      tSSq[t - 1] += dat**2

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for t in range(Npts - 1):
    Kdat[t].append(tK[t] / float(count))
    KdatSq[t].append(tKSq[t] / float(count))
    Sdat[t].append(tS[t] / float(25. * count))      # Average over all
    SdatSq[t].append(tSSq[t] / float(25. * count))  # Average over all
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

times = np.array(range(Nt - 1))
Ktot = np.array([sum(Kdat[t]) for t in range(Npts - 1)])
KtotSq = np.array([sum(KdatSq[t]) for t in range(Npts - 1)])
Stot = np.array([sum(Sdat[t]) for t in range(Npts - 1)])
StotSq = np.array([sum(SdatSq[t]) for t in range(Npts - 1)])

# Monitor how many times the fits don't converge
NK = 0;     Kfails = 0
NS = 0;     Sfails = 0

# Fit results for all jk samples
Kout = np.empty((len(K_in), Nblocks), dtype = np.float)
Sout = np.empty((len(S_in), Nblocks), dtype = np.float)
for i in range(Nblocks):    # Jackknife samples
  K = np.empty(Nt - 1, dtype = np.float)
  Kerr = np.empty(Nt - 1, dtype = np.float)
  S = np.empty(Nt - 1, dtype = np.float)
  Serr = np.empty(Nt - 1, dtype = np.float)
  for t in range(Npts - 1):
    K[t] = (Ktot[t] - Kdat[t][i]) / (Nblocks - 1.)
    temp = (KtotSq[t] - KdatSq[t][i]) / (Nblocks - 1.)
    Kerr[t] = np.sqrt(temp - K[t]**2)
    S[t] = (Stot[t] - Sdat[t][i]) / (Nblocks - 1.)
    temp = (StotSq[t] - SdatSq[t][i]) / (Nblocks - 1.)
    Serr[t] = np.sqrt(temp - S[t]**2)
  # Account for t <--> Nt - t symmetry
  for t in range(Npts - 1, Nt - 1):
    K[t] = K[Nt - 1 - t];
    Kerr[t] = Kerr[Nt - 1 - t];
    S[t] = S[Nt - 1 - t];
    Serr[t] = Serr[Nt - 1 - t];

  # Fits!
  temp, cov, infodict, mesg, success = optimize.leastsq(Kerrfunc, K_in[:], args=(times, K, Kerr), full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    for j in range(len(K_in)):
      Kout[j][i] = np.fabs(temp[j])   # Cosh is symmetric under -m<-->m
    NK += 1
  else:                               # Ignore failed fits
    print mesg
    Kfails += 1

  temp, cov, infodict, mesg, success = optimize.leastsq(Serrfunc, S_in[:], args=(times, S, Serr), full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    for j in range(len(S_in)):
      Sout[j][i] = np.fabs(temp[j])   # Cosh is symmetric under -m<-->m
    NS += 1
  else:                               # Ignore failed fits
    print mesg
    Sfails += 1
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
# Ignore any failed fits
print "Fitting with %d blocks of length %d MDTU..." % (Nblocks, block_size)
print "%d of %d Konishi fits succeeded" % (NK, Nblocks)
print "%d of %d SUGRA fits succeeded" % (NS, Nblocks)

outfile = open('results/fit.corr', 'w')
print >> outfile, "# Fitting with %d blocks of length %d MDTU" \
                  % (Nblocks, block_size)
print >> outfile, "%d of %d Konishi fits succeeded" % (NK, Nblocks)
print >> outfile, "%d of %d SUGRA fits succeeded" % (NS, Nblocks)

if NK > 0:
  tag = ["constant", "amplitude", "mass"]
  for i in range(len(K_in)):
    ave = np.average(Kout[i])
    var = (NK - 1.) * np.sum((Kout[i] - ave)**2) / float(NK)
    print >> outfile, "Konishi %s %.6g %.4g" % (tag[i], ave, np.sqrt(var))

if NS > 0:
  tag = ["amplitude", "mass"]
  for i in range(len(S_in)):
    ave = np.average(Sout[i])
    var = (NS - 1.) * np.sum((Sout[i] - ave)**2) / float(NS)
    print >> outfile, "SUGRA %s %.6g %.4g" % (tag[i], ave, np.sqrt(var))

runtime += time.time()
print >> outfile, "# Runtime: %.2g seconds" % runtime

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
