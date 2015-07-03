#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Monitor progress and performance of pfaffian measurements

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

# Check progress and performance of each measurement
for i in range(len(toCheck)):
  complete = -1
  done = 0
  matvec_per_sec = 0.0
  count = 0
  for line in open(toCheck[i]):
    # Format: Q has # columns --> # matvecs and # MBytes per core...
    if line.startswith('Q has '):
      total = int((line.split())[5])

    # Format: Columns #--# of #: # matvecs in # seconds re im
    elif "matvecs in" in line:
      temp = line.split()
      done += int(temp[4])
      matvec_per_sec += float(temp[4]) / float(temp[7])
      count += 1

    elif line.startswith('RUNNING COMPLETED'):
      matvec_per_sec /= count     # Average matvecs per second
      print toCheck[i], "is already complete",
      print "(performed %.1f matvecs per second)" % matvec_per_sec
      complete = 1

  if complete < 0:
    matvec_per_sec /= count     # Average matvecs per second
    todo = total - done         # Matvecs still to do
    time = todo / (float)(matvec_per_sec * 3600.0)
    print toCheck[i], "should finish in %.1f hours" % time,
    print "(performing %.1f matvecs per second," % matvec_per_sec,
    percent = 100.0 * done / float(total)
    print "%.0f%% through)" % percent
# ------------------------------------------------------------------
