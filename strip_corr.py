#!/usr/bin/python
import os
import sys
import glob
from subprocess import call
# ------------------------------------------------------------------
# Strip cheap and in some cases incorrect results out of existing
# (unsmeared) measurement files
# Currently strips the following list of tags:
tags = ['KONISHI ', 'SUGRA ', 'd_correlator_r: ', 'CORR_', 'rsymm: ', \
        'INVLINK ', 'RSYMM ', 'hvy_pot', 'PLOT_LOOP ', 'DL_LOOP ', \
        'PLOLAR_LOOP ', 'Fixing ', 'GFIX ', 'BEFORE ', 'AFTER ', \
        'POT_LOOP ', 'D_LOOP ', 'POLAR_LOOP ', 'FLAVOR ', 'Error ']

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Cycle over measurement files
for filename in glob.glob('Out/corr.*'):
  outfile = open('TEMP', 'w')
  check = -1
  for line in open(filename):
    toPrint = 1
    for tag in tags:
      if line.startswith(tag):
        toPrint = -1
        break   # No need to check rest of tags
    if toPrint > 0:
      print >> outfile, line.rstrip()
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print filename, "did not complete"
    sys.exit(1)
  outfile.close()
  os.rename('TEMP', filename)
# ------------------------------------------------------------------

