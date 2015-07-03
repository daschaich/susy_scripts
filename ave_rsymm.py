#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Print blocked plaquette and R-symmetry-modified Wilson loop averages
# with blocked standard errors

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
# Third argument tells us whether to analyze "corr" or "stout" files
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
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
  if cfg not in cfgs and cfg > cut:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files", files, "found"
  sys.exit(1)

# If we're missing some initial measurements,
# increase thermalization cut
cut = cfgs[0]

# Determine largest Wilson loops from first output file
firstfile = 'Out/' + tag + '.' + str(cfgs[0])
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)

MAX = -1
for line in open(firstfile):
  if line.startswith('N=4 SYM,'):
    Nc = int((((line.split())[4]).split(','))[0])
  elif line.startswith('rsymm: MAX = '):
    MAX = int((line.split())[3])
  elif line.startswith('INVLINK '):
    break   # Done scanning through file

if MAX < 0:
  print "ERROR:", firstfile, "did not use up-to-date R symmetry measurement"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct arrays of blocked measurements for each observable
# Consider 1x1, ..., 1xMAX, 2x1, ..., MAXx1, ..., MAXxMAX
# The second index (plus 1) gives the number of inverted links
# W = usual Wilson loop, M = modified loop
# D = difference (M - W), R = ratio (M - W)/(M + W)
Wdat = [[[] for i in range(MAX)] for i in range(MAX)]
Mdat = [[[] for i in range(MAX)] for i in range(MAX)]
Ddat = [[[] for i in range(MAX)] for i in range(MAX)]
Rdat = [[[] for i in range(MAX)] for i in range(MAX)]

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0         # How many measurements in each block
begin = cut       # Where each block begins, to be incremented

# Accumulators
tW = np.zeros((MAX, MAX), dtype = np.float)
tM = np.zeros((MAX, MAX), dtype = np.float)
tD = np.zeros((MAX, MAX), dtype = np.float)
tR = np.zeros((MAX, MAX), dtype = np.float)
for MDTU in cfgs:
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    if count == 0:
      print "ERROR: no data to average after file %s:" % toOpen
      sys.exit(1)
    for i in range(MAX):
      for j in range(MAX):
        Wdat[i][j].append(tW[i][j] / float(20. * count))
        Mdat[i][j].append(tM[i][j] / float(20. * count))
        Ddat[i][j].append(tD[i][j] / float(20. * count))
        Rdat[i][j].append(tR[i][j] / float(20. * count))
        tW[i][j] = 0.0
        tM[i][j] = 0.0
        tD[i][j] = 0.0
        tR[i][j] = 0.0
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)

  # Running averages
  filename = 'Out/' + tag + '.' + str(MDTU)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen
  check = -1
  polar = -1    # !!! Temporary check for Konishi, which should probably be reversed for R symmetry transformations
  for line in open(toOpen[0]):
    # Skip measurement if average link looks unreasonably large
    # This seems to be due to occasional near-singular link matrix inversions
#    if line.startswith('INVLINK '):
#      temp = float((line.split())[6])
#      if temp > Nc:
#        print "WARNING: skipping %s due to large link inverse %.4g" \
#              % (toOpen[0], temp)
#        check = 1   # Avoid spurious non-completion error
#        break

    # Format: RSYMM normal [dir] inverted [dir] usual mod
    # This comes last in the output files,
    # so it will use the correct sav
    if line.startswith('RSYMM '):
      temp = line.split()
      norm = int(temp[1]) - 1
      inv = int(temp[3]) - 1
      W_dat = float(temp[5])
      M_dat = float(temp[6])
      D_dat = (M_dat - W_dat) / float(2.0 * inv + 2.0)
      tW[norm][inv] += W_dat
      tM[norm][inv] += M_dat
      tD[norm][inv] += D_dat
      tR[norm][inv] += D_dat / (M_dat + W_dat)
      # Only tick counter once per measurement
      if norm == 0 and inv == 0 and temp[2] == '[0]' and temp[4] == '[1]':
        count += 1

    elif line.startswith('rsymm considering polar-projected links'):
      polar = 1

    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print infile, "reports two measurements"
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)
  if polar == -1:
    print toOpen[0], "did not consider polar-projected links"
    sys.exit(1)

# Check special case that last block is full
# Assume last few measurements are equally spaced
if cfgs[-1] >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  for i in range(MAX):
    for j in range(MAX):
      if count == 0:
        print "ERROR: no data to average after file %s:" % toOpen
        sys.exit(1)
      Wdat[i][j].append(tW[i][j] / float(20. * count))
      Mdat[i][j].append(tM[i][j] / float(20. * count))
      Ddat[i][j].append(tD[i][j] / float(20. * count))
      Rdat[i][j].append(tR[i][j] / float(20. * count))
  # Record block data
  block_data[0].append(count)
  block_data[1].append(begin)
  block_data[2].append(begin + block_size)

Nblocks = len(Wdat[0][0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
if Nblocks == 1:
  print "ERROR: need multiple blocks to take average"
  sys.exit(1)

print "Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
outfile = open('results/rsymm.dat', 'w')
print >> outfile, "# Averaging with %d blocks of length %d MDTU" % (Nblocks, block_size)
print >> outfile, "# norm inv diff err rel err modified err usual err"

for i in range(MAX):
  for j in range(MAX):
    print >> outfile, "%d %d" % (i + 1, j + 1),

    dat = np.array(Ddat[i][j])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
    print >> outfile, "%.6g %.4g" % (ave, err),

    dat = np.array(Rdat[i][j])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
    print >> outfile, "%.6g %.4g" % (ave, err),

    dat = np.array(Mdat[i][j])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
    print >> outfile, "%.6g %.4g" % (ave, err),

    dat = np.array(Wdat[i][j])
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(Nblocks - 1.)
    print >> outfile, "%.6g %.4g" % (ave, err)

# More detailed block information
#for i in range(Nblocks):
#  print >> outfile, \
#        "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
outfile.close()
# ------------------------------------------------------------------
