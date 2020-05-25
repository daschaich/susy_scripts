#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import emcee.autocorr as acor
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# with given thermalization cut and block size

# Assume one ensemble per directory (overwrite results files)
# Assume Polyakov loop data are already normalized by Nc

# Parse arguments: first is thermalization cut, second is block size
# Will check block size is larger than poly_mod, 9th scalar square
# and lowest eigenvalue auto-correlation times
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print("Usage:", str(sys.argv[0]), "<cut> <block>")
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print("ERROR: data/ does not exist")
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
if t_cut < (0.7 * float(cut)):
  t_block /= 2

final_MDTU = float(temp[1])
if good == -1:
  print("Error: no data to analyze ", end='')
  print("since cut=%d but we only have %d MDTU" % (cut, final_MDTU))
  sys.exit(1)

# We are good to go.  Save the path, which may be useful to flag errors
path = os.getcwd()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Check that block size is larger than poly_mod, 9th scalar square
# and lowest eigenvalue auto-correlation times
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
    print("Error: poly_mod meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep))
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[1]))

# Arguments explained in emcee.readthedocs.io/en/stable/user/autocorr/
#                    and emcee.readthedocs.io/en/stable/tutorials/autocorr/
# 'c' is step size for window search (default 5)
#     Larger c should reduce noise, but can add bias...
# 'tol' is min ratio between data and tau (default 50)
# 'Quiet' prints warning rather than shutting down if tol not satisfied
tau = acor.integrated_time(np.array(dat), c=5, tol=10, quiet=True)
tau *= sep
if tau > block_size:
  print("Error: poly_mod autocorrelation time %d " % tau, end='')
  print("is larger than block size %d " % block_size, end='')
  print("in %s" % path)
  sys.exit(1)

# Record poly_mod auto-correlation time for future reference
# Include effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/poly_mod.autocorr'
outfile = open(outfilename, 'w')
print("%d # %d" % (tau, eff_stat), file=outfile)
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
    print("Error: scalarsquares meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep))
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[-1]))

# Arguments discussed above
tau = acor.integrated_time(np.array(dat), c=5, tol=10, quiet=True)
tau *= sep
if tau > block_size:
  print("Error: scalar square autocorrelation time %d " % tau, end='')
  print("is larger than block size %d " % block_size, end='')
  print("in %s" % path)
  sys.exit(1)

# Record scalar square auto-correlation time for future reference
# Include effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/scalarsquares.autocorr'
outfile = open(outfilename, 'w')
print("%d # %d" % (tau, eff_stat), file=outfile)
outfile.close()

# Finally, for the lowest eigenvalue we have only one datum per line
# Format: MDTU,min_eig
dat = []
sep = 10      # !!!Now non-trivial
prev = 0
for line in open('data/eig.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = int(temp[0])

  # Need to check separation and update prev before skipping to therm cut
  if not MDTU - prev == sep:
    print("Error: min_eig meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep))
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[-1]))

# Arguments discussed above
tau = acor.integrated_time(np.array(dat), c=5, tol=10, quiet=True)
tau *= sep
if tau > block_size:
  print("Error: lowest eigenvalue autocorrelation time %d " % tau, end='')
  print("is larger than block size %d " % block_size, end='')
  print("in %s" % path)
  sys.exit(1)

# Record lowest eigenvalue auto-correlation time for future reference
# Include effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/eig.autocorr'
outfile = open(outfilename, 'w')
print("%d # %d" % (tau, eff_stat), file=outfile)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the Polyakov loop, bosonic action, fermion action,
# 'Myers' scalar trilinear term and SO(6) / SO(3) action ratio
# we're interested in the first datum on each line (following MDTU)
# For the Polyakov loop, this is the (Nc-normalized) modulus
for obs in ['poly_mod', 'SB', 'SF', 'Myers', 'ratio']:
  ave = 0.0         # Accumulate within each block
  count = 0
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

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      ave += float(temp[1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        datList.append(ave / float(count))
        begin += block_size
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For algorithmic/cost quantities
# we're interested in the first datum on each line
# Need to work in terms of trajectories rather than MDTU
for obs in ['wallTU', 'cg_iters', 'accP', 'exp_dS']:
  ave = 0.0         # Accumulate within each block
  count = 0
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

    # Accumulate within block
    elif traj > begin and traj <= (begin + t_block):
      ave += float(temp[1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if traj == (begin + t_block):
        datList.append(ave / float(count))
        begin += t_block
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif traj > (begin + t_block):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
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
    print("WARNING: No scalar eigenvalue data")
    continue

  ave = []      # Accumulate within each block
  count = 0
  datList = []
  for i in range(Nc):
    ave.append(0.0)
    datList.append([])

  begin = cut       # Where each block begins, to be incremented
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue

    # Accumulate within each block
    elif MDTU > begin and MDTU <= (begin + block_size):
      for i in range(Nc):
        ave[i] += float(temp[i + 1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        for i in range(len(ave)):
          datList[i].append(ave[i] / float(count))
          ave[i] = 0.0
        begin += block_size
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(Nc):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave[i] = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print("%.8g %.4g " % (ave[i], err), end='', file=outfile)
  print("# %d" % N, file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the scalar squares we're interested in 9 data on each line
# Also monitor 6-->4+2 splitting due to discretization artifacts
for obs in ['scalarsquares']:
  ave = []          # Accumulate within each block
  count = 0
  datList = []

  # For splitting
  high = 0.0
  low  = 0.0
  mid  = 0.0
  split = 0.0
  splitList = []
  Nscalar = 9
  for i in range(Nscalar):
    ave.append(0.0)
    datList.append([])

  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue

    # Accumulate within each block
    elif MDTU > begin and MDTU <= (begin + block_size):
      high = float(temp[8]) + float(temp[9])
      low  = float(temp[4]) + float(temp[5]) + float(temp[6]) + float(temp[7])
      mid  = (high + low) / 6.0
      for i in range(Nscalar):
        ave[i] += float(temp[i + 1])
      split += (0.5 * high - 0.25 * low) / mid
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        splitList.append(split / float(count))
        split = 0.0
        for i in range(len(ave)):
          datList[i].append(ave[i] / float(count))
          ave[i] = 0.0
        begin += block_size
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(Nscalar):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave[i] = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print("%.8g %.4g " % (ave[i], err), end='', file=outfile)
  print("# %d" % N, file=outfile)
  outfile.close()

  # Print splitting separately
  dat = np.array(splitList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/splitting.dat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------
