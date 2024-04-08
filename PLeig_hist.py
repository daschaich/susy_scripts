#!/usr/bin/env python3
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Polyakov loop eigenvalues
# Use Out/eig files that include these measurements
# Superimpose Ndirs=3 ensembles specified by input arguments
# Assume only thermalized measurements have been run
# Save resulting plot as ./PLeig_hist.pdf

# Parse arguments: Directories for SU(8), SU(12) and SU(16) ensembles,
# the upper bound for the y-axis,
# and finally a tag for the title for the plot
if len(sys.argv) < 6:
  print("Usage:", str(sys.argv[0]), end='')
  print("<SU(8) dir> <SU(12) dir> <SU(16) dir>", end='')
  print("<y-axis upper bound> <plot title tag>")
  sys.exit(1)

Ndirs = 3
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
#  tocheck = dirname + '/Out/PLeig.*'
  tocheck = dirname + '/Out/eig.*'
  files = glob.glob(tocheck)
  if len(files) == 0:
    print("ERROR: no files", tocheck)
    sys.exit(1)

  # Extract Nc from dirname -- it's after 'Nc' then before '_'
  temp = (dirname.split('Nc'))[1]
  Nc = int((temp.split('_'))[0])

  # Cycle through output files to load data
  for filename in files:
    check = -1
    for line in open(filename):
      # Format: LINES_EIG {Nc x phase}
      if line.startswith('LINES_EIG '):
        temp = line.split()
        for i in range(Nc):
          dat[count].append(float(temp[-1 - i]))
      elif line.startswith('RUNNING COMPLETED'):
        check = 1
    if check == -1:
      print(filename, "did not complete")
      sys.exit(1)

  # Check that we have all the data we should
  if not len(dat[count]) == Nc * len(files):
    print("ERROR: Have %d data from %d SU(%d) measurements" \
          % (len(dat[count]), len(files), int(Nc)))
    sys.exit(1)

  # Check that all phases are within [-pi, pi),
  # accounting for rounding in output files
  # Can be moved into parse_BMN.py (cf. parse_SYM.py)
  if max(dat[count]) > 3.142 or min(dat[count]) < -3.142:
    print("ERROR: %s phases exceed [-pi, pi): %.4g %.4g" \
          % (dirname, max(dat[count]), min(dat[count])))
    sys.exit(1)
  count += 1

# Create histogram
nbins = 10
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat[0], nbins, log=False, density=True, align='mid',
         edgecolor='red', label='SU(8)', histtype='step', hatch='||')
plt.hist(dat[1], nbins, log=False, density=True, align='mid',
         edgecolor='green', label='SU(12)', histtype='step', hatch='\\\\')
plt.hist(dat[2], nbins, log=False, density=True, align='mid',
         edgecolor='blue', label='SU(16)', histtype='step', hatch='//')

plt.axis([-np.pi, np.pi, 0.0, ymax])
plt.xticks([-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi],
           ['$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', \
                      r'$\frac{\pi}{2}$', r'$\pi$'])
plt.grid(False)

plt.title('Phase of Polyakov loop eigenvalues, ' + title)
plt.xlabel('Phase')
plt.ylabel('Relative frequency')
plt.legend()

# Save a pdf
outfile = 'PLeig_hist.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
