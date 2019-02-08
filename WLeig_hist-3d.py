#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Wilson line eigenvalues from output files 'WLeig'
# Superimpose multiple ensembles specified by input arguments
# For now only consider unitarized (and not full) Wilson lines in z-dir
# Assume only thermalized measurements have been run
# Save resulting plot as ./WLeig_hist.pdf

# Parse arguments: L=8 and 12 ensembles,
# the upper bound for the y-axis,
# and finally a tag for the title for the plot
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<L=8 dir> <L=12 dir>",
  print "<y-axis upper bound> <plot title tag>",
  sys.exit(1)

Ndirs = 2
dirnames = []
for i in range(Ndirs):
  dirnames.append(str(sys.argv[i + 1]))
ymax = float(sys.argv[-2])
title = str(sys.argv[-1])

# Prepare data to plot from each directory
count = 0
dat = [[] for x in range(Ndirs)]
for dirname in dirnames:
  # Check that we have output files to analyze
  tocheck = dirname + '/Out/WLeig.*'
  files = glob.glob(tocheck)
  if len(files) == 0:
    print "ERROR: no files", tocheck
    sys.exit(1)

  # Extract Nc from dirname -- it's after 'Nc' then before '_'
  temp = (dirname.split('Nc'))[1]
  Nc = int((temp.split('_'))[0])

  # Cycle through output files to load data
  for filename in files:
    check = -1
    for line in open(filename):
      if line.startswith('nt '):
        nt = float((line.split())[1])
      # Format: LINES_POLAR_EIG x y z t dir {Nc x phase}
      elif line.startswith('LINES_POLAR_EIG '):
        temp = line.split()
        for i in range(Nc):
          dat[count].append(float(temp[-1 - i]))
      elif line.startswith('RUNNING COMPLETED'):
        check = 1
    if check == -1:
      print filename, "did not complete"
      sys.exit(1)

  # Check that we have all the data we should
  if not len(dat[count]) == Nc * nt * len(files):
    print "ERROR: Have %d data from %d SU(%d) measurements with L=%d" \
          % (len(dat[count]), len(files), int(Nc), nt)
    sys.exit(1)

  # Check that all phases are within [-pi, pi),
  # accounting for rounding in output files
  # Can be commented out to speed up analysis
  if max(dat[count]) > 3.142 or min(dat[count]) < -3.142:
    print "ERROR: %s phases exceed [-pi, pi): %.4g %.4g" \
          % (dirname, max(dat[count]), min(dat[count]))
    sys.exit(1)
  count += 1

# Create histogram
nbins = 10
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat[0], nbins, log=False, normed=True, align='mid',
         edgecolor='blue', label='L=8', histtype='step', hatch='//')
plt.hist(dat[1], nbins, log=False, normed=True, align='mid',
         edgecolor='green', label='L=12', histtype='step', hatch='\\\\')

plt.axis([-np.pi, np.pi, 0.0, ymax])
plt.xticks([-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi],
           ['$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', \
                      r'$\frac{\pi}{2}$', r'$\pi$'])
plt.grid(False)

plt.title('Phase of unitarized Wilson line eigenvalues, ' + title)
plt.xlabel('Phase')
plt.ylabel('Relative frequency')
plt.legend()

# Save a pdf
outfile = 'WLeig_hist.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
