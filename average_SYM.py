#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# given a thermalization cut and block size
# Assume one ensemble per directory
# Assume Polyakov loop data are properly normalized by Nc
# ------------------------------------------------------------------



# ------------------------------------------------------------------
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
MDTUfile = 'data/TU.csv'
good = -1
for line in open(MDTUfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[1]) > cut:
    good = 1
    break

final_MDTU = float(temp[1])
if good == -1:
  print "Error: no data to analyze",
  print "since cut=%d but we only have %d MDTU" % (cut, final_MDTU)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette is special -- average two data per line
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
plaqfile = 'data/plaq.csv'
for line in open(plaqfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU < cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave += (float(temp[1]) + float(temp[2])) / 2
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next block
    datList.append(ave / count)
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = (float(temp[1]) + float(temp[2])) / 2

# Now print mean and standard error, assuming N>1
dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
outfilename = 'results/plaq.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the Polyakov loop, bosonic action, average link and determinant
# we're interested in the first datum on each line
# For the Polyakov loop, this is the (Nc-normalized) modulus
# For the determinant, this is |det - 1|
for obs in ['poly_mod', 'SB', 'Flink', 'det']:
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU < cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif float(temp[0]) >= (begin + block_size):  # Move on to next block
      datList.append(ave / count)
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])

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
