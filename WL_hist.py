#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Wilson line magnitudes from output files 'WLeig'
# For now only consider unitarized (and not full) Wilson lines
# Assume only thermalized measurements have been run
# Save resulting plot as ./WL_hist_$tag.pdf
# Extract $tag from path---no input arguments

if len(sys.argv) > 1:
  print "Usage:", str(sys.argv[0])
  sys.exit(1)

# Check that we have output files to analyze
files = glob.glob('Out/WLeig.*')
if len(files) == 0:
  print "ERROR: no files Out/WLeig.* in this directory"
  sys.exit(1)

# Extract Nc from path -- it's after 'Nc' then before '_'
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = int((temp.split('_'))[0])

# Extract tag from path as everything after the last '/'
tag = (cwd.split('/'))[-1]

# Cycle through output files to load data---normalize by Nc
dat = []
for filename in files:
  check = -1
  for line in open(filename):
    if line.startswith('nt '):
      nt = float((line.split())[1])
    # Format: LINES_POLAR_TR  x y z t dir real imag
    elif line.startswith('LINES_POLAR_TR  '):
      temp = line.split()
      re = float(temp[-2])
      im = float(temp[-1])
      dat.append(np.sqrt(re * re + im * im) / float(Nc))
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print filename, "did not complete"
    sys.exit(1)

# Check that we have all the data we should
if not len(dat) == nt * len(files):
  print "ERROR: Have %d data from %d SU(%d) measurements with L=%d" \
        % (len(dat), len(files), int(Nc), nt)
  sys.exit(1)

# Create histogram
nbins = 20
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, bins=np.arange(0.0, 1.0, 0.05),
         log=False, normed=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, 1.0])
plt.grid(False)

plt.xlabel('|WL|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
outfile = 'WL_hist_' + tag + '.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
