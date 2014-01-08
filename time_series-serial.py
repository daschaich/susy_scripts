#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Average over SUGRA correlator components to extract time series data
# with proper normalization

# Parse argument: the time component to extract
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<time>"
  sys.exit(1)
time = int(sys.argv[1])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile('corrSUGRA'):
  print "ERROR: corrSUGRA does not exist"
  sys.exit(1)

# Average over SUGRA components for each measurement
component = 0
toOpen = 'corrSUGRA'
toAve = 0.
for line in open(toOpen):
  # If we're done with a measurement, print it and reset for the next
  if component == 25:
    toAve /= 20.  # Five diagonal components already subtracted
    print toAve
    toAve = 0.
    component = 0

  # Running average
  temp = line.split()
  if int(temp[0]) == time:
    toAve += float(temp[1])
    component += 1
# ------------------------------------------------------------------
