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
# Check that block size is larger than poly_mod auto-correlation time
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
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette is special -- average two data per line
# Also look at difference plaq_ss - plaq_st
ave = 0.0         # Accumulate within each block
diff = 0.0
count = 0
datList = []
diffList = []
begin = cut       # Where each block begins, to be incremented
plaqfile = 'data/plaq.csv'
for line in open(plaqfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue

  # Accumulate within each block
  elif MDTU > begin and MDTU <= (begin + block_size):
    ave += 0.5 * (float(temp[1]) + float(temp[2]))
    diff += float(temp[1]) - float(temp[2])
    count += 1

    # If that "<=" is really "==" then we are done with this block
    # Record it and re-initialize for the next block
    if MDTU == (begin + block_size):
      datList.append(ave / float(count))
      diffList.append(diff / float(count))

      begin += block_size
      ave = 0.0
      diff = 0.0
      count = 0

  # This doesn't happen for ensembles I generate
  # May need to be revisited for more general applicability
  elif MDTU > (begin + block_size):
    print "ERROR: Unexpected behavior in %s, aborting" % obsfile
    sys.exit(1)

# Now print mean and standard error, assuming N>1
dat = np.array(datList, dtype = np.float64)
N = np.size(dat)
ave = np.mean(dat)
err = np.std(dat) / np.sqrt(N - 1.0)
outfilename = 'results/plaq.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()

dat = np.array(diffList, dtype = np.float64)
ave = np.mean(dat)
err = np.std(dat) / np.sqrt(N - 1.0)
outfilename = 'results/plaq_diff.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the Polyakov loop, bosonic action, fermion action, average link,
# and monopole world line density
# we're interested in the first datum on each line
# For the Polyakov loop, this is the (Nc-normalized) modulus
for obs in ['poly_mod', 'poly_mod_polar', 'SB', 'SF', 'Flink', 'mono']:
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

    # Accumulate within each block
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
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
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

    # Accumulate within each block
    elif traj > begin and traj <= (begin + t_block):
      ave += float(temp[1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if traj >= (begin + t_block):
        datList.append(ave / float(count))

        begin += t_block
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif traj > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the plaquette determinant we're interested in all three data on each line
# These are |det-1|^2, 1-Re(det) and Im(det)
for obs in ['det']:
  ave = [0.0, 0.0, 0.0]       # Accumulate within each block
  count = 0
  datList = [[], [], []]
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue

    elif MDTU > begin and MDTU <= (begin + block_size):
      ave[0] += float(temp[1])
      ave[1] += float(temp[2])
      ave[2] += float(temp[3])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        for i in range(len(ave)):
          datList[i].append(ave[i] / float(count))

        begin += block_size
        for i in range(len(ave)):
          ave[i] = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(len(ave)):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave[i] = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print >> outfile, "%.8g %.4g" % (ave[i], err),
  print >> outfile, "# %d" % N
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the widths and Wilson lines
# we're interested in all four data on each line
# These widths are for plaq, Re(det), Im(det) and Tr[U.Udag]/N
for obs in ['widths', 'lines_mod', 'lines_mod_polar']:
  ave = [0.0, 0.0, 0.0, 0.0]      # Accumulate within each block
  count = 0
  datList = [[], [], [], []]
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
      ave[0] += float(temp[1])
      ave[1] += float(temp[2])
      ave[2] += float(temp[3])
      ave[3] += float(temp[4])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        for i in range(len(ave)):
          datList[i].append(ave[i] / float(count))

        begin += block_size
        for i in range(len(ave)):
          ave[i] = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(len(ave)):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave[i] = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print >> outfile, "%.8g %.4g" % (ave[i], err),
  print >> outfile, "# %d" % N
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

        begin += block_size
        for i in range(len(ave)):
          ave[i] = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  for i in range(Nc):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave[i] = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print >> outfile, "%.8g %.4g" % (ave[i], err),
  print >> outfile, "# %d" % N
  outfile.close()
# ------------------------------------------------------------------

