#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compute pfaffian from printed diagonal elements of Q, as a sanity check

# Parse argument: the file(s) to check
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toCheck = []
for i in range(1, len(sys.argv)):
  toCheck.append(str(sys.argv[i]))

# First make sure we're calling this from the right place
for i in range(len(toCheck)):
  if not os.path.isfile(toCheck[i]):
    print "ERROR: %s does not exist" % toCheck[i]
    sys.exit(1)

# Average over SUGRA components for each measurement
phase = 0.0
log_mag = 0.0
for i in range(len(toCheck)):
  check = -1
  for line in open(toCheck[i]):
    # Format: Columns #--# of #: pivoting #<--># to obtain re im
    #         Columns #--# of #: no need to pivot to obtain re im
    if line.startswith('Columns '):
      temp = line.split()
      re = float(temp[-2])
      im = float(temp[-1])

      log_mag += np.log(re**2 + im**2) / 2.0
      phase += np.angle(re + 1j * im)
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toCheck[i], "did not complete"
    sys.exit(1)

tr = np.fmod(phase, 2.0 * np.pi)
if tr < 0:
  phase = tr + 2.0 * np.pi
else:
  phase = tr

print "%.8g %.8g %.8g %.8g" % \
      (log_mag, phase, np.fabs(np.cos(phase)), np.fabs(np.sin(phase)))
# ------------------------------------------------------------------
