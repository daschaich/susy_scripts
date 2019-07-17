#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data file to compute susceptibility
# for the unitarized Polyakov loop
# (Other targets may be added in the future)

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

# Extract Nc from path -- it's after 'Nc' then before '_'
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = int((temp.split('_'))[0])
norm = Nc * Nc

final_MDTU = float(temp[1])
if good == -1:
  print "Error: no data to analyze",
  print "since cut=%d but we only have %d MDTU" % (cut, final_MDTU)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# We're interested in the first datum on each line
# This is the modulus (followed by the real and imaginary parts)
for obs in ['poly_mod']:
  skip = -1
  count = 0
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  datList = []
  sqList = []
  begin = cut       # Where each block begins, to be incremented

  obsfile = 'data/' + obs + '.csv'
  if not os.path.isfile(obsfile):
    print "ERROR: " + obsfile + " does not exist"
    sys.exit(1)

  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue
    elif MDTU > begin and MDTU < (begin + block_size):
      tr = float(temp[1])
      ave += tr
      aveSq += tr * tr
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d MDTU" % (obs, int(MDTU))
        skip = 1
        break
      datList.append(ave / count)
      sqList.append(aveSq / count)
      begin += block_size
      count = 1                         # Next block begins here
      tr = float(temp[1])
      ave = tr
      aveSq = tr * tr

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print "ERROR: need multiple blocks to take average"
    sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  #   chi = N^2 * [<PL^2> - <PL>^2]
  dat = np.array(datList)
  tot = sum(dat)
  sq = np.array(sqList)
  totSq = sum(sq)
  chi = np.zeros(N, dtype = np.float)
  for i in range(N):    # Jackknife samples
    vev = (tot - dat[i]) / (N - 1.0)
    sq_vev = (totSq - sq[i]) / (N - 1.0)
    chi[i] = sq_vev - vev * vev

  # Now we can average over jackknife samples and print out results
  ave = np.mean(chi, dtype = np.float64)
  var = (N - 1.0) * np.sum((chi - ave)**2) / float(N)
  outfilename = 'results/' + obs + '.suscept'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (norm * ave, norm * np.sqrt(var), N)
  outfile.close()
# ------------------------------------------------------------------

