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
matvecs = 0    # Make sure all columns are accounted for
phase = 0.0
log_mag = 0.0
total = -1    # In case we're just checking a subset without "Q has # columns"
for i in range(len(toCheck)):
  check = -1
  for line in open(toCheck[i]):
    # Format: Q has # columns --> # matvecs and # MBytes per core...
    if line.startswith('Q has '):
      temp = line.split()
      total = int(temp[5])

    # Format: Columns #--# of #: # matvecs in # seconds re im
    elif line.startswith('Columns '):
      temp = line.split()
      matvecs += int(temp[4])
      re = float(temp[9])
      im = float(temp[10])

      log_mag += np.log(re**2 + im**2) / 2.0
      phase += np.angle(re + 1j * im)
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toCheck[i], "did not complete"
    sys.exit(1)

if total > 0 and matvecs != total:
  print "ERROR: Only counted %d of %d matvecs" % (matvecs, total)

tr = np.fmod(-1.0 * phase, 2.0 * np.pi)
if tr < 0:
  phase = tr + 2.0 * np.pi
else:
  phase = tr

print "%.8g %.8g %.8g %.8g" % \
      (-1.0 * log_mag, phase, np.fabs(np.cos(phase)), np.fabs(np.sin(phase)))
# ------------------------------------------------------------------
