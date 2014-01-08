#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy.misc import comb   # For N choose k
# ------------------------------------------------------------------
# Compute cumulants of Wilson loop data

# Parse arguments: first is thermalization cut,
# remainder is tag to match, such as "POT_LOOP 1 1 1 " or "D_LOOP   1 1 1 "
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> '<tag>'"
  sys.exit(1)
cut = int(sys.argv[1])
tag = str(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct list of which configurations have been analyzed
cfgs = []
allfiles = glob.glob('Out/corr.*')
for filename in allfiles:
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

# Extract maximum temporal extent of loop from first file
for line in open(allfiles[0]):
  if line.startswith('hvy_pot'):
    temp = line.split()
    # Format: hvy_pot(): MAX_T = #, MAX_X = #
    Nt = int((temp[3].split(','))[0])   # Cut comma from end
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Load all log(data) into arrays
# Shift data from [1, Nt] to [0, Nt - 1]
# Discard complete measurements if any of the data are negative
datarray = [[] for x in range(Nt)]

for MDTU in cfgs:
  filename = 'Out/corr.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "# WARNING: multiple files named %s:" % filename,
    print toOpen
  for line in open(toOpen[0]):
    # Format: LABEL x y z t dat
    if line.startswith(tag):
      temp = line.split()
      t = int(temp[4]) - 1            # Shift Nt
      dat = float(temp[5])
      if dat < 0:
        print "# WARNING: skipping %s due to negative datum %.4g for t=%d" \
              % (filename, dat, t + 1)
        for i in range(t):    # Clear already-recorded data
          datarray[i].pop()
        break                 # Moves us on to next MDTU
      datarray[t].append(np.log(dat))

# Check that all the arrays are the same length
Ndat = len(datarray[0])
for t in range(1, Nt):
  if len(datarray[t]) != Ndat:
    print "ERROR: data mismatch for t=%d: %d should be %d" \
          % (t, len(datarray[t]), Ndat)
    sys.exit(1)
print "# Analyzing %d of %d total measurements" % (Ndat, len(cfgs))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Jackknife the first few moments and cumulants
if Ndat == 1:
  print "ERROR: not enough data for jackknifing"
  sys.exit(1)

Norder = 5   # How many cumulants to consider

# All JK results
muJK = np.empty((Norder, Nt, Ndat), dtype = np.float)
kaJK = np.empty((Norder, Nt, Ndat), dtype = np.float)
for t in range(Nt):
  alldat = np.array(datarray[t])    # All data for this t
  for i in range(Ndat):
    dat = np.delete(alldat, i)    # Remove ith data point
    for n in range(Norder):
      muJK[n][t][i] = np.mean(dat**(n + 1), dtype = np.float64)
      kaJK[n][t][i] = muJK[n][t][i]
      for m in range(n):
        # comb(N, k) is N-choose-k
        # Have checked that (np.std(dat))**2 reproduces ka[2]
        kaJK[n][t][i] -= comb(n, m) * kaJK[m][t][i] * muJK[n - m - 1][t][i]

  # Average over jackknife samples for ka_2 through ka_N
  print t + 1,
  for n in range(1, Norder):
    ave = np.average(kaJK[n][t])
    var = (Ndat - 1.) * np.sum((kaJK[n][t] - ave)**2) / float(Ndat - 1.)
    print "%.4g %.4g" % (ave, np.sqrt(var)),
  print
# ------------------------------------------------------------------
