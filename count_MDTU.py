#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Count MDTU in "*out*" files in the current directory
# Print the running total

# Construct list of out* files
cfgs = []
for filename in glob.glob('*out*'):
  cfg = int(filename.split('.')[-1])    # Number after last .
  if cfg not in cfgs:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: no files *out* found"
  sys.exit(1)

# Accumulate running some from ordered list of files
totMDTU = 0
for i in cfgs:
  filename = '*out*.' + str(i)
  toOpen = glob.glob(filename)
  if len(toOpen) > 1:
    print "ERROR: multiple files named %s:" % filename,
    print toOpen

  # We have a file, let's cycle through it's lines
  done = -1
  for line in open(toOpen[0]):
    if line.startswith('traj_length '):
      length = float((line.split())[1])
    elif line.startswith('traj_between_meas '):
      N = int((line.split())[1])
    if line.startswith('RUNNING COMPLETED'):
      done = 1    # This job actually finished

  if done > 0:
    totMDTU += N * length
    print "%g MDTU after file %s" % (totMDTU, toOpen[0])
  else:
    print toOpen[0], "did not complete"
    sys.exit(1)
# ------------------------------------------------------------------
