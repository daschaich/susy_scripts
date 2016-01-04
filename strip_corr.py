#!/usr/bin/python
import os
import sys
import glob
from subprocess import call
# ------------------------------------------------------------------
# Strip cheap and in some cases incorrect results out of existing
# (unsmeared) measurement files
# Currently strips the following list of tags:
tags = ['KONISHI ', 'SUGRA ', 'd_correlator_r: ', 'CORR_', 'VACSUB_', 'vevK' \
        'rsymm: ', 'INVLINK ', 'RSYMM ', \
        'hvy_pot', 'PLOT_LOOP ', 'DL_LOOP ', 'PLOLAR_LOOP ', \
        'Fixing ', 'GFIX ', 'BEFORE ', 'AFTER ', \
        'POT_LOOP ', 'D_LOOP ', 'POLAR_LOOP ', 'FLAVOR ', 'Error ', \
        'LINES ', 'LINES_POLAR ', 'UUBAR_EIG ', 'POLAR_EIG ']

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Cycle over measurement files
for filename in glob.glob('Out/corr.*'):
  outfile = open('TEMP', 'w')
  check = -1
  for line in open(filename):
    if line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print filename, "reports two measurements"
      check = 1

    toPrint = 1
    for tag in tags:
      if line.startswith(tag):
        toPrint = -1
        break   # No need to check rest of tags
    if toPrint > 0:
      print >> outfile, line.rstrip()

  if check == -1:
    print filename, "did not complete"
    continue
  outfile.close()
  os.rename('TEMP', filename)
# ------------------------------------------------------------------

