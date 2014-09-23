#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Blocked Konishi or SUGRA correlator C(t) and Delta C / Delta t
# projected to zero spatial momentum

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument is the file to analyze
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <file>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
toOpen = str(sys.argv[3])

# Make sure given file exists
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Figure out Nt from given file
times = []
for line in open(toOpen):
  t = int(line.split()[0])
  if len(times) > 0 and t == times[0]:
    break
  else:
    times.append(t)

# We can ignore all t > Nt / 2
Nt = len(times) / 2 + 1
if (times[-1] + 1) / 2 + 1 != Nt:
  print "Something funny going on:", times
  sys.exit(1)

# Now construct arrays of blocked measurements
corr = [[] for x in range(Nt)]
diff = [[] for x in range(Nt - 1)]
meas = 1
prev_time = float('NaN')    # To make it obvious if this isn't overwritten
begin = cut       # Where each block begins, to be incremented
toAve = [0 for x in range(Nt)]
dAve = [0 for x in range(Nt - 1)]
for line in open(toOpen):
  # If we're done with a block, record it and reset for the next
  if meas >= (begin + block_size):
    for t in range(Nt):
      corr[t].append(toAve[t] / float(block_size))
    for t in range(Nt - 1):
      diff[t].append(dAve[t] / float(block_size))
    # Record and reset block data
    begin += block_size
    toAve = [0 for x in range(Nt)]
    dAve = [0 for x in range(Nt - 1)]

  # Running averages
  temp = line.split()
  time = int(temp[0])
  if time >= Nt:
    continue
  if meas >= begin and meas < (begin + block_size):
    dat = float(temp[1])
    toAve[time] += dat
    if time > 0:
      dAve[time - 1] += dat - prev_time
    prev_time = dat

  # Increment meas after entire measurement recorded
  if time == (Nt - 1):
    meas += 1

Nblocks = len(corr[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "# Averaging with %d blocks of length %d measurements" \
      % (Nblocks, block_size)
print "# C(t):"
for t in range(Nt):
  dat = np.array(corr[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print "%d %.6g %.4g" % (t, ave, err)

print "# Delta C / Delta t:"
for t in range(Nt - 1):
  dat = np.array(diff[t])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
  print "# %.2g %.6g %.4g" % (float(t) + 0.5, ave, err)
# ------------------------------------------------------------------
