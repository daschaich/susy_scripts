#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data file to compute specific heat
# for the bosonic BMN model
# !!! As with bBMN_package.py, fixing Nt=24

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
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

# Extract Nc from path --- it's after 'Nc' then before '_'
Nt = 24.0
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = float((temp.split('_'))[0])

# Extract t from path --- it's everything after 't'
t = float((cwd.split('t'))[-1])
cube_root = 1.0 / (Nt * t)              # lambda_lat^{1/3}

# Overall normalization factor
norm = cube_root * cube_root * Nt * Nt / (Nc * Nc)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# We're interested in both data on each line
# The first is E/N^2 while the second is E_prime
# We want to compute <E^2> - <E>^2 - <E_prime>, with no 1/N^2 factors
for obs in ['energy']:
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  aveP = 0.0        # For E_prime
  count = 0
  datList = []
  sqList = []
  pList = []
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
      tr = Nc * Nc * float(temp[1])     # Convert from E/N^2
      ave += tr
      aveSq += tr * tr
      aveP += float(temp[2])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU >= (begin + block_size):
        datList.append(ave / count)
        sqList.append(aveSq / count)
        pList.append(aveP / count)
        begin += block_size
        ave = 0.0
        aveSq = 0.0
        aveP = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif traj > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  #   C_v = (lalat^{2/3} Nt^2 / Nc^2) [<E^2> - <E>^2 - <E_prime>]
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  tot = sum(dat)
  sq = np.array(sqList, dtype = np.float64)
  totSq = sum(sq)
  p = np.array(pList, dtype = np.float64)
  totP = sum(p)
  Cv = np.zeros_like(dat)
  for i in range(N):    # Jackknife samples
    vev = (tot - dat[i]) / (N - 1.0)
    sq_vev = (totSq - sq[i]) / (N - 1.0)
    p_vev = (totP - p[i]) / (N - 1.0)
    Cv[i] = sq_vev - vev * vev - p_vev

  # Now we can average over jackknife samples and print out results
  ave = np.mean(Cv)
  var = (N - 1.0) * np.sum((Cv - ave)**2) / float(N)
  outfilename = 'results/' + obs + '.specheat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (norm * ave, norm * np.sqrt(var), N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------

