#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi and SUGRA correlator averages
# with blocked standard errors

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
if not os.path.isfile('corrK'):
  print "ERROR: corrK does not exist"
  sys.exit(1)
if not os.path.isfile('corrSUGRA'):
  print "ERROR: corrSUGRA does not exist"
  sys.exit(1)

# Extract lattice volume from path
# For now only Nt is used; we assume it is a two-digit number!!!
path = os.getcwd()
temp = path.split('/')
Nt = int(temp[-2][-2:])     # Last two digits before "/r#"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# K = Konishi, S = SUGRA
Kdat = [[] for x in range(Nt)]
Sdat = [[] for x in range(Nt)]

# Konishi is fairly straightforward
meas = 1
begin = cut       # Where each block begins, to be incremented
toOpen = 'corrK'
toAve = [0 for x in range(Nt)]
for line in open(toOpen):
  # If we're done with a block, record it and reset for the next
  if meas >= (begin + block_size):
    for t in range(Nt):
      # Correct normalization issue
      Kdat[t].append(toAve[t] / (25. * float(block_size)))
    # Record and reset block data
    begin += block_size
    toAve = [0 for x in range(Nt)]

  # Running averages
  temp = line.split()
  time = int(temp[0])
  if meas >= begin and meas < (begin + block_size):
    toAve[time] += float(temp[1])

  # Increment meas after entire measurement recorded
  if time == (Nt - 1):
    meas += 1

# SUGRA has complication of 25 components for each measurement
meas = 1
component = 1
begin = cut       # Where each block begins, to be incremented
toOpen = 'corrSUGRA'
toAve = [0 for x in range(Nt)]
for line in open(toOpen):
  # If we're done with a block, record it and reset for the next
  if meas >= (begin + block_size):
    for t in range(Nt):
      Sdat[t].append(toAve[t] / (25. * float(block_size)))
    # Record and reset block data
    begin += block_size
    toAve = [0 for x in range(Nt)]

  # Running averages
  temp = line.split()
  time = int(temp[0])
  if meas >= begin and meas < (begin + block_size):
    toAve[time] += float(temp[1])

  # Increment component after every Nt lines
  if time == (Nt - 1):
    component += 1
    # Increment meas after every 25 components
    if component == 26:
      meas += 1
      component = 1

Nblocks = len(Kdat[0])
if len(Sdat[0]) != Nblocks:
  print "ERROR: different counts for Konishi (%d) and SUGRA (%d)..." \
        % (Nblocks, len(Sdat[0]))
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d measurements" \
      % (Nblocks, block_size)
print "Konishi:"
for t in range(Nt):
  dat = np.array(Kdat[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print "%d %.6g %.4g" % (t, ave, err)

print "\nSUGRA:"
for t in range(Nt):
  dat = np.array(Sdat[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print "%d %.6g %.4g" % (t, ave, err)
# ------------------------------------------------------------------
