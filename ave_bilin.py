#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked fermion bilinear and Ward identity violation averages
# with blocked standard errors

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us whether to analyze "corr" or "stout" files
if len(sys.argv) < 4:
  print("Usage:", str(sys.argv[0]), "<cut> <block> <tag>")
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print("ERROR: Out/ does not exist")
  sys.exit(1)

# Try to set C2 from path
C2 = 1.0
temp = os.getcwd()
if '-c' in temp:
  C2 = float((temp.split('-c'))[1])   # Everything after '-c'

# Construct list of which configurations have been analyzed
cfgs = []
files = 'Out/' + tag + '.*'
for filename in glob.glob(files):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs and cfg > cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print("ERROR: no files", files, "found")
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each observable
# Ignore stochastic noise in measurements
# Need to add factor of C2 to gauge piece
# F = fermion bilinear piece, G = gauge piece,
# D = susy transform (difference C2 * G - F),
# N = normalized (C2 * G - F) / (C2 * G + F)
Fdat = []
Gdat = []
Ddat = []
Ndat = []

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tF = 0.0
tG = 0.0
tD = 0.0
tN = 0.0
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    Fdat.append(tF / float(count))
    Gdat.append(tG / float(count))
    Ddat.append(tD / float(count))
    Ndat.append(tN / float(count))
    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    tF = 0.0
    tG = 0.0
    tD = 0.0
    tN = 0.0

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print("ERROR: multiple files named %s:" % filename, end='', file=outfile)
    print(toOpen)
  check = -1
  for line in open(toOpen[0]):
    # Format: SUSY f_dat imag g_dat diff [with C2=1]
    # Imaginary component should average to zero
    if line.startswith('SUSY '):
      count += 1
      temp = line.split()
      diff = C2 * float(temp[3]) - float(temp[1])
      tF += float(temp[1])
      tG += C2 * float(temp[3])
      tD += diff
      a = float(temp[1])
      b = C2 * float(temp[3])
      tN += diff / np.sqrt(a * a + b * b)
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print(toOpen[0], "did not complete")
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  Fdat.append(tF / float(count))
  Gdat.append(tG / float(count))
  Ddat.append(tD / float(count))
  Ndat.append(tN / float(count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Fdat)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print("ERROR: need multiple blocks to take average")
  sys.exit(1)

print("Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size))
outfile = open('results/bilin.dat', 'w')
print("# Averaging with %d blocks of length %d MDTU" \
      % (Nblocks, block_size), file=outfile)
print("# diff err rel err gauge err bilin err", file=outfile)

for obs in [Ddat, Ndat, Gdat, Fdat]:
  dat = np.array(obs)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.0)
  print("%.6g %.4g " % (ave, err), end='', file=outfile)
print("# %d" % Nblocks, file=outfile)

# More detailed block information
#for i in range(Nblocks):
#  print("# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i]), \
#        file=outfile)
outfile.close()
# ------------------------------------------------------------------
