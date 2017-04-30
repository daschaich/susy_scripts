#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Wilson line eigenvalues from output files 'WLeig'
# For now only consider unitarized (and not full) Wilson lines
# Assume only thermalized measurements have been run
# Save resulting plot as WLeig_hist_$tag.pdf

if len(sys.argv) > 1:
  print "Usage:", str(sys.argv[0])
  sys.exit(1)

# Check that we have output files to analyze
dat = []
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

# Cycle through output files to load data
for filename in files:
  check = -1
  for line in open(filename):
    if line.startswith('nt '):
      nt = float((line.split())[1])
    # Format: LINES_POLAR_EIG x y t dir Nc*phase
    elif line.startswith('LINES_POLAR_EIG '):
      temp = line.split()
      for i in range(Nc):
        dat.append(float(temp[-1 - i]))
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print filename, "did not complete"
    sys.exit(1)

# Check that we have all the data we should
if not len(dat) == Nc * nt * len(files):
  print "ERROR: Have %d data from %d SU(%d) measurements with L=%d" \
        % (len(dat), len(files), int(Nc), nt)
  sys.exit(1)

# Create histogram
nbins = 20
n, bins, patches = plt.hist(dat, nbins, log=False, normed=True, align='mid',
                            facecolor='blue', alpha=0.9, label=tag,
                            histtype='stepfilled')

plt.axis([-np.pi, np.pi, 0.0, 1.0])
plt.xticks([-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi],
           ['$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', \
                      r'$\frac{\pi}{2}$', r'$\pi$'])
plt.grid(True)

plt.title('Phase of unitarized Wilson line eigenvalues')
plt.xlabel('Phase')
plt.ylabel('Relative frequency')
plt.legend()

# Save a pdf
outfile = 'WLeig_hist_' + tag + '.pdf'
plt.savefig(outfile)
# ------------------------------------------------------------------
