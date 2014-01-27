#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compute pfaffian from printed diagonal elements of Q, as a sanity check

# Parse argument: the file to check
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toCheck = str(sys.argv[1])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toCheck):
  print "ERROR: %s does not exist" % toCheck
  sys.exit(1)

# Average over SUGRA components for each measurement
phase = 0.
log_mag = 0.
for line in open(toCheck):
  if line.startswith('Columns '):
    temp = line.split()
    re = float(temp[9])
    im = float(temp[10])

    log_mag += np.log(re**2 + im**2) / 2.0
    phase += np.angle(re + 1j * im)

tr = np.fmod(-1.0 * phase, 2.0 * np.pi)
if tr < 0:
  phase = tr + 2.0 * np.pi
else:
  phase = tr

print "%.8g %.8g %.8g %.8g" % \
      (-1.0 * log_mag, phase, np.fabs(np.cos(phase)), np.fabs(np.sin(phase)))
# ------------------------------------------------------------------
