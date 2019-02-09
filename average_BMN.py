#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import acor         # Uses "the Kubo formula" to compute autocorrelation time
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# given a thermalization cut and block size
# Assume one ensemble per directory
# Assume Polyakov loop data are properly normalized by Nc

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Check that we actually have data to average
# and convert thermalization cut from MDTU to trajectory number
MDTUfile = 'data/TU.csv'
sav = 0
good = -1
for line in open(MDTUfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[1]) > cut:
    good = 1
    t_cut = sav
    break
  sav = float(temp[0])

# Guess whether we should also convert the block size
# from MDTU to trajectory number
# They differ when tau=2 trajectories are used...
t_block = block_size
if t_cut < float(cut) / 1.5:
  t_block /= 2

final_MDTU = float(temp[1])
if good == -1:
  print "Error: no data to analyze",
  print "since cut=%d but we only have %d MDTU" % (cut, final_MDTU)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Check that block size is larger
# than poly_mod and 9th scalar square auto-correlation times
# TODO: Look at smallest eigenvalue as well?
# For poly_mod we want the first datum on each line (following MDTU)
# Format: MDTU,|Tr(L)|,ReTr(L),ImTr(L)
dat = []
sep = 1
prev = 0
for line in open('data/poly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = int(temp[0])

  # Need to check separation and update prev before skipping to therm cut
  if not MDTU - prev == sep:
    print "Error: poly_mod meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep)
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[1]))

# Discard this mean and sigma -- we'll recompute it later
tau, mean, sigma = acor.acor(np.array(dat))
tau *= sep
if tau > block_size:
  print "Error: poly_mod autocorrelation time %d" % tau,
  print "is larger than block size %d" % block_size,
  print "in %s" % path
  sys.exit(1)

# Record poly_mod auto-correlation time for future reference
# Include average and effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/poly_mod.autocorr'
outfile = open(outfilename, 'w')
print >> outfile, "%d --> %.8g %.4g # %d" % (tau, mean, sigma, eff_stat)
outfile.close()

# Next, for the scalar square we want the last (tenth) datum on each line
# Format: MDTU,Tr(X1)^2,...,Tr(X9)^2
dat = []
sep = 1       # Redundant, retained to be explicit
prev = 0
for line in open('data/scalarsquares.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = int(temp[0])

  # Need to check separation and update prev before skipping to therm cut
  if not MDTU - prev == sep:
    print "Error: scalarsquares meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep)
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[-1]))

# Discard this mean and sigma -- we'll recompute it later
tau, mean, sigma = acor.acor(np.array(dat))
tau *= sep
if tau > block_size:
  print "Error: scalar square autocorrelation time %d" % tau,
  print "is larger than block size %d" % block_size,
  print "in %s" % path
  sys.exit(1)

# Record scalar square auto-correlation time for future reference
# Include average and effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/scalarsquares.autocorr'
outfile = open(outfilename, 'w')
print >> outfile, "%d --> %.8g %.4g # %d" % (tau, mean, sigma, eff_stat)
outfile.close()
# ------------------------------------------------------------------




# ------------------------------------------------------------------
# For the Polyakov loop, bosonic action, fermion action,
# energy and 'Myers' scalar trilinear term
# we're interested in the first datum on each line (following MDTU)
# For the Polyakov loop, this is the (Nc-normalized) modulus
for obs in ['poly_mod', 'SB', 'SF', 'energy', 'Myers']:
  skip = -1
  count = 0
  ave = 0.0         # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue
    elif MDTU > begin and MDTU < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d MDTU" % (obs, int(MDTU))
        skip = 1
        break
      datList.append(ave / count)
      begin += block_size
      count = 1                         # Next block begins here
      ave = float(temp[1])

  if len(datList) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For algorithmic/cost quantities
# we're again interested in the first datum on each line
# but have to work in terms of trajectories rather than MDTU
for obs in ['wallTU', 'cg_iters', 'accP', 'exp_dS']:
  skip = -1
  count = 0
  ave = 0.0         # Accumulate within each block
  datList = []
  begin = t_cut     # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M') or line.startswith('t'):
      continue
    temp = line.split(',')
    traj = float(temp[0])
    if traj <= t_cut:
      continue
    elif traj > begin and traj < (begin + t_block):
      ave += float(temp[1])
      count += 1
    elif traj >= (begin + t_block):     # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d traj" % (obs, int(traj))
        skip = 1
        break
      datList.append(ave / count)
      begin += t_block
      count = 1                         # Next block begins here
      ave = float(temp[1])

  if len(datList) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# For the scalar eigenvalues we're interested in Nc data on each line
for obs in ['scalar_eig_ave']:
  # Figure out Nc from number of points on first non-trivial line
  Nc = -1
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    Nc = len(temp) - 1
    break

  if Nc < 0:
    print "WARNING: No scalar eigenvalue data"
    continue

  ave = []      # Accumulate within each block
  datList = []
  for i in range(Nc):
    ave.append(0.0)
    datList.append([])

  skip = -1
  count = 0
  begin = cut       # Where each block begins, to be incremented
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue
    elif MDTU > begin and MDTU < (begin + block_size):
      for i in range(Nc):
        ave[i] += float(temp[i + 1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d MDTU" % (obs, int(MDTU))
        skip = 1
        break
      for i in range(len(ave)):
        datList[i].append(ave[i] / count)
        ave[i] = float(temp[i + 1])     # Next block begins here
      begin += block_size
      count = 1

  if len(datList[0]) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(Nc):
    dat = np.array(datList[i])
    N = np.size(dat)
    ave[i] = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
    print >> outfile, "%.8g %.4g" % (ave[i], err),
  print >> outfile, "# %d" % N
  outfile.close()
# ------------------------------------------------------------------

