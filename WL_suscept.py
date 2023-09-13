#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data file to compute susceptibility
# for the unitarized spatial Wilson line
# (Other targets may be added in the future)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 4:
  print("Usage:", str(sys.argv[0]), "<cut> <block> <dir>")
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
direction = int(sys.argv[3]) # x, y, z for 1, 2, 3, respectively
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print("ERROR: data/ does not exist")
  sys.exit(1)

# Extract Nc from path -- it's after 'Nc' then before '_'
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = int((temp.split('_'))[0])
norm = Nc * Nc
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# We're interested in the given datum on each line
# For dir=1,2,3 this is the x,y,z-direction modulus, respectively
# For dir=5 we combine the y- and z-directions
for obs in ['lines_mod_polar']:
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  count = 0
  datList = []
  sqList = []
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
      if direction < 4:
        tr = float(temp[direction])
        ave += tr
        aveSq += tr * tr
        count += 1
      elif direction == 5:
        tr = float(temp[2])
        ave += tr
        aveSq += tr * tr
        tr = float(temp[3])
        ave += tr
        aveSq += tr * tr
        count += 2
      else:
        print("Error: Unrecognized dir %d" % direction)

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU >= (begin + block_size):
        datList.append(ave / count)
        sqList.append(aveSq / count)
        begin += block_size
        ave = 0.0
        aveSq = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif traj > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  #   chi = N^2 * [<WL^2> - <WL>^2]
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  tot = sum(dat)
  sq = np.array(sqList, dtype = np.float64)
  totSq = sum(sq)
  chi = np.zeros_like(dat)
  for i in range(N):    # Jackknife samples
    vev = (tot - dat[i]) / (N - 1.0)
    sq_vev = (totSq - sq[i]) / (N - 1.0)
    chi[i] = sq_vev - vev * vev

  # Now we can average over jackknife samples and print out results
  ave = np.mean(chi)
  var = (N - 1.0) * np.sum((chi - ave)**2) / float(N)
  outfilename = 'results/' + obs + '.suscept'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (norm * ave, norm * np.sqrt(var), N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------

