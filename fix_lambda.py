#!/usr/bin/python
import os
import sys
import glob
import math
from subprocess import call
# ------------------------------------------------------------------
# Impose consistent conventions for lambda and error_per_site,
# which were temporarily redefined for "slnc_code_2" computations
# While lambda is safe to recompute arbitrarily many times,
# we need to make sure that error_per_site is only replaced once

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Extract Nc from path -- it's after 'Nc' then before '_'
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = float((temp.split('_'))[0])

# Cycle over all output files (including eig, phase, etc)
for filename in glob.glob('Out/*'):
  outfile = open('TEMP', 'w')
  check = -1
  for line in open(filename):
    # Check for slnc_code vs. slnc_code_2
    if 'slnc_code' in line:
      if 'slnc_code_2' in line:   # Good to go
        check = 1
      else:                       # Move to next file
        break

    if line.startswith('lambda '):
      continue
    # Need to make sure that error_per_site is only replaced once
    # Fortunately it always seems to have been 1e-5
    elif line.startswith('error_per_site 1e-05'):
      err = math.sqrt(float((line.split())[1]))
      print >> outfile, "error_per_site %g" % err
    elif line.startswith('lambda='):
      temp = line.rstrip()    # Extract kappa from end of line
      ka = float((temp.split('='))[-1])
      la = 0.5 * Nc / ka
      print >> outfile, "kappa=Nc/(2lambda)=%g --> lambda=%g" % (ka, la)
    else:
      print >> outfile, line.rstrip()

  outfile.close()
  if check > 0:
    os.rename('TEMP', filename)
  else:
    os.remove('TEMP', filename)
# ------------------------------------------------------------------

