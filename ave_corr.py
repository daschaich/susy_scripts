#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked Konishi, SUGRA and Konishi superpartner correlator averages
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
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct list of which configurations have been analyzed
cfgs = []
for filename in glob.glob('Out/corr.*'):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs and cfg > cut:
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

# We exploit t <--> Nt - t symmetry
# To only print 0 through Nt / 2 (a total of Npts points)
Npts = Nt / 2 + 1  # Assume Nt is even
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each correlator
# For now just average over all 25 SUGRA components
# K = Konishi, S = SUGRA
Kdat = [[] for x in range(Npts)]
Sdat = [[] for x in range(Npts)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tK = [0 for x in range(Npts)]
tS = [0 for x in range(Npts)]
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    for t in range(Npts):
      Kdat[t].append(tK[t] / float(count))
      Sdat[t].append(tS[t] / float(25. * count))  # Average over all
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tK = [0 for x in range(Npts)]
    tS = [0 for x in range(Npts)]

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
      time = int(temp[1])
      if time == 0:
        count += 1    # Only increment once per measurement!
      dat = float(temp[2])
      tK[time] += dat
    # Format: SUGRA a b t dat
    elif line.startswith('SUGRA '):
      temp = line.split()
      time = int(temp[3])
      dat = float(temp[4])
      tS[time] += dat

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for t in range(Npts):
    Kdat[t].append(tK[t] / float(count))
    Sdat[t].append(tS[t] / float(25. * count))  # Average over all
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Kdat[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

aveK = [0 for x in range(Nt)]
errK = [0 for x in range(Nt)]
aveS = [0 for x in range(Nt)]
errS = [0 for x in range(Nt)]
for t in range(Npts):
  # Konishi
  dat = np.array(Kdat[t])
  aveK[t] = np.mean(dat, dtype = np.float64)
  errK[t] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)

  # SUGRA
  dat = np.array(Sdat[t])
  aveS[t] = np.mean(dat, dtype = np.float64)
  errS[t] = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)

# Account for t <--> Nt - t symmetry
for t in range(Npts, Nt):
  aveK[t] = aveK[Nt - t];
  errK[t] = errK[Nt - t];
  aveS[t] = aveS[Nt - t];
  errS[t] = errS[Nt - t];

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Kfile = open('results/konishi.dat', 'w')
print >> Kfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
Sfile = open('results/sugra.dat', 'w')
print >> Sfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
for t in range(Nt):
  print >> Kfile, "%d %.6g %.4g # %d" % (t, aveK[t], errK[t], Nblocks)
  print >> Sfile, "%d %.6g %.4g # %d" % (t, aveS[t], errS[t], Nblocks)

# More detailed block information
#for i in range(Nblocks):
#  for outfile in [Kfile, Sfile]:
#    print >> outfile, \
#          "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#          % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
Kfile.close()
Sfile.close()
# ------------------------------------------------------------------
