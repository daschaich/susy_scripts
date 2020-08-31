#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Wilson line eigenvalues
# Use data/WLeig.csv file for time-series plots
# Superimpose five ensembles specified by input arguments
#   First three (filled) are N=4, 6 and 8 with L=12
#   Last two (unfilled) are L=8 and 16 with N=8
# For now only consider unitarized (and not full) Wilson lines in z-dir
# !!! Assume only thermalized measurements have been run
# Save resulting plot as ./WLeig_hist.pdf

# Parse arguments: L=8, 12 and 16 ensembles,
# the upper bound for the y-axis,
# and finally a tag for the title for the plot
if len(sys.argv) < 8:
  print("Usage:", str(sys.argv[0]), end='')
  print("<U(4) dir> <U(6) dir> <U(8) dir>", end='')
  print("<L=8 dir> <L=16 dir>", end='')
  print("<y-axis upper bound> <plot title tag>")
  sys.exit(1)

Ndirs = 5
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
  infile = dirname + '/data/WLeig.csv'
  if not os.path.isfile(infile):
    print("Problem opening", infile)
    continue    # Skip to next file
    sys.exit(1)

  # Extract Nc from dirname -- it's after 'Nc' then before '_'
  temp = (dirname.split('Nc'))[1]
  Nc = int((temp.split('_'))[0])

  # Go through time-series data files to load data
  # Nc phases on each line --- parse_SYM.py checks all within [-pi, pi)
  for line in open(infile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    for i in range(Nc):
      dat[count].append(float(temp[1 + i]))
  count += 1

# Create histogram
nbins = 10
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat[0], nbins, log=False, density=True, align='mid',
         edgecolor='red', label='N=4, L=12', histtype='step', hatch='||')
plt.hist(dat[1], nbins, log=False, density=True, align='mid',
         edgecolor='green', label='N=6, L=12', histtype='step', hatch='\\\\')
plt.hist(dat[2], nbins, log=False, density=True, align='mid',
         edgecolor='blue', label='N=8, L=12', histtype='step', hatch='//')
plt.hist(dat[3], nbins, log=False, density=True, align='mid',
         edgecolor='orange', label='N=8, L=8', histtype='step')
plt.hist(dat[4], nbins, log=False, density=True, align='mid',
         edgecolor='purple', label='N=8, L=16', histtype='step')

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
