#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Run jackknifed fits of zero-momentum-projected Konishi and SUGRA correlators
# C_K(t) = C + A(e^{-M_K t} + e^{-M_K (Nt - t)})
# C_S(t) = A(e^{-M_S t} + e^{-M_S (Nt - t)})

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

# Now we can define the fit functions for C_K(t) and C_S(t)
# errfunc will be minimized via least-squares optimization
Kcorr = lambda p, t: p[0] + p[1] * (np.exp(-p[2] * t) + np.exp(-p[2] * (Nt - 1 - t)))
Kerrfunc = lambda p, t, y, err: (Kcorr(p, t) - y) / err
K_in = [20., 0.0001, 1.]          # Order-of-magnitude initial guesses
Scorr = lambda p, t: p[0] * (np.exp(-p[1] * t) + np.exp(-p[1] * (Nt - 1 - t)))
Serrfunc = lambda p, t, y, err: (Scorr(p, t) - y) / err
S_in = [0.001, 0.1]               # Order-of-magnitude initial guesses

# We exploit t <--> Nt - t symmetry to only print 0 through Nt / 2
Npts = Nt / 2 + 1 - tmin  # Assume Nt is even
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# There is some funny business with indexing from zero despite tmin>0
# We need to shift things by tmin to deal with it
# For now just average over all 25 SUGRA components
# K = Konishi, S = SUGRA
Kdat   = [[] for x in range(Npts)]
KdatSq = [[] for x in range(Npts)]
Sdat   = [[] for x in range(Npts)]
SdatSq = [[] for x in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK   = [0 for x in range(Npts)]
tKSq = [0 for x in range(Npts)]
tS   = [0 for x in range(Npts)]
tSSq = [0 for x in range(Npts)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for t in range(Npts):
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
    tK   = [0 for x in range(Npts)]
    tKSq = [0 for x in range(Npts)]
    tS   = [0 for x in range(Npts)]
    tSSq = [0 for x in range(Npts)]

  # Running averages
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
      if t < tmin:
        continue
      elif t == tmin:
        count += 1    # Only increment once per measurement!
      dat = float(temp[2])
      tK[t - tmin] += dat
      tKSq[t - tmin] += dat**2
    # Format: SUGRA a b t dat
    elif line.startswith('SUGRA '):
      temp = line.split()
      t = int(temp[3])
      if t < tmin:
        continue
      dat = float(temp[4])
      tS[t - tmin] += dat
      tSSq[t - tmin] += dat**2

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for t in range(Npts):
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

times  = np.arange(tmin, Nt / 2 + 1, dtype=np.int)
Ktot   = np.array([sum(Kdat[t]) for t in range(Npts)])
KtotSq = np.array([sum(KdatSq[t]) for t in range(Npts)])
Stot   = np.array([sum(Sdat[t]) for t in range(Npts)])
StotSq = np.array([sum(SdatSq[t]) for t in range(Npts)])

# Monitor how many times the fits don't converge
NK = 0;     Kfails = 0
NS = 0;     Sfails = 0

# Fit results for all jk samples
#np.seterr(all='raise')
Kout = np.empty((len(K_in), Nblocks), dtype = np.float)
Sout = np.empty((len(S_in), Nblocks), dtype = np.float)
for i in range(Nblocks):    # Jackknife samples
  K    = np.empty(Npts, dtype = np.float)
  Kerr = np.empty(Npts, dtype = np.float)
  S    = np.empty(Npts, dtype = np.float)
  Serr = np.empty(Npts, dtype = np.float)
  for t in range(Npts):
    K[t] = (Ktot[t] - Kdat[t][i]) / (Nblocks - 1.)
    temp = (KtotSq[t] - KdatSq[t][i]) / (Nblocks - 1.)
    Kerr[t] = np.sqrt(temp - K[t]**2)
    S[t] = (Stot[t] - Sdat[t][i]) / (Nblocks - 1.)
    temp = (StotSq[t] - SdatSq[t][i]) / (Nblocks - 1.)
    Serr[t] = np.sqrt(temp - S[t]**2)

  # Fits!
#  success = -1    # In case fit fails at runtime
#  temp, cov, infodict, mesg, success = optimize.leastsq(Kerrfunc, K_in[:],
#                                                        args=(times, K, Kerr),
#                                                        full_output=True)
#  if success in (1, 2, 3, 4):         # Fit succeeded
#    for j in range(len(K_in)):
#      Kout[j][i] = temp[j]
#    NK += 1
#  else:                               # Ignore failed fits
#    print mesg
#    Kfails += 1

  success = -1    # In case fit fails at runtime
  temp, cov, infodict, mesg, success = optimize.leastsq(Serrfunc, S_in[:],
                                                        args=(times, S, Serr),
                                                        full_output=True)
  if success in (1, 2, 3, 4):         # Fit succeeded
    for j in range(len(S_in)):
      Sout[j][i] = temp[j]
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
