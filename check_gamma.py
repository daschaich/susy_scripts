#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Compute gamma = mu / sqrt(lambda) = mu Nt / rt
# and (mu N_t)^2 / lambda = (N_t gamma)^2

# Extract rt from path -- it's after 'rt' then before '_'
cwd = os.getcwd()
temp = (cwd.split('rt'))[1]
rt = float((temp.split('_'))[0])

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Grab last output file
infile = (sorted(glob.glob('Out/out.*')))[-1]
for line in open(infile):
  if line.startswith('nz '):
    L = float((line.split())[-1])
  elif line.startswith('nt '):
    Nt = float((line.split())[-1])
  elif line.startswith('lambda '):
    la = float((line.split())[-1])
  elif line.startswith('bmass '):
    mu = float((line.split())[-1])
  elif line.startswith('kappa='):
    temp = line.rstrip()    # Extract lambda from end of line
    la = float((temp.split('='))[-1])

rx = rt * L / Nt
#print "rt, rx = %g, %g" % (rt, rx)
#print "lambda = %g" % la
#print "mu = %g" % mu

ga = mu / math.sqrt(la)
check = mu * Nt / rt
#print "gamma = %.4g = %.4g" % (ga, check)

combo = (mu * Nt)**2 / la
check = (Nt * ga)**2
#print "scale = %.4g = %.4g" % (combo, check)

# Format for tex tables...
print "%.3f & %.3f & $%d\X %d$  & %.7f & %.7f & %.2f                                              & %.2f" \
      % (rt, rx, int(L), int(Nt), la, mu, ga, combo)
# ------------------------------------------------------------------

