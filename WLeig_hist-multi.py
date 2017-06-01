#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Wilson line eigenvalues from output files 'WLeig'
# Superimpose multiple ensembles specified by input arguments
# For now only consider unitarized (and not full) Wilson lines
# Assume only thermalized measurements have been run
# Save resulting plot as ./WLeig_hist.pdf

# Parse arguments: first is number of ensembles,
# each of which much then be listed
# The last argument is the title for the plot
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<# of dirs>",
  print "<name> of each dir",
  sys.exit(1)
Ndirs = int(sys.argv[1])
dirnames = []
for i in range(Ndirs):
  dirnames.append(str(sys.argv[i + 1]))
title = str(sys.argv[-1])

# Prepare data to plot from each directory
count = 0
dat = [[] for x in range(Ndirs)]
for dirname in dirnames:
  # Check that we have output files to analyze
  files = glob.glob(dirname + 'Out/WLeig.*')
  if len(files) == 0:
    print "ERROR: no files Out/WLeig.* in this directory"
    sys.exit(1)

  # Extract Nc from path -- it's after 'Nc' then before '_'
  cwd = os.getcwd()
  temp = (cwd.split('Nc'))[1]
  Nc = int((temp.split('_'))[0])

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

# Create histogram
nbins = 15
plt.hist(dat[0], nbins, log=False, normed=True, align='mid',
         facecolor='blue', alpha=0.9, label='SU(6)', histtype='stepfilled')
plt.hist(dat[1], nbins, log=False, normed=True, align='mid',
         facecolor='green', alpha=0.9, label='SU(9)', histtype='stepfilled')
plt.hist(dat[2], nbins, log=False, normed=True, align='mid',
         facecolor='red', alpha=0.9, label='SU(12)', histtype='stepfilled')

plt.axis([-np.pi, np.pi, 0.0, 1.0])
plt.xticks([-np.pi, -0.5 * np.pi, 0.0, 0.5 * np.pi, np.pi],
           ['$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', \
                      r'$\frac{\pi}{2}$', r'$\pi$'])
plt.grid(True)

plt.title('Phase of unitarized Wilson line eigenvalues,' + title)
plt.xlabel('Phase')
plt.ylabel('Relative frequency')
plt.legend()

# Save a pdf
outfile = 'WLeig_hist.pdf'
plt.savefig(outfile)
# ------------------------------------------------------------------
