#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Polyakov loop magnitudes
# Save resulting plot as ./PL_hist_$tag.pdf
# Extract $tag from path rather than input argument

# Parse argument: Thermalization cut
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<cut>"
  sys.exit(1)
cut = int(sys.argv[1])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Extract tag from path as everything after the last '/'
cwd = os.getcwd()
tag = (cwd.split('/'))[-1]

# Extract Nc-normalized modulus as third datum on each line of data file
# Format: MDTU,|Tr(L)|,ReTr(L),ImTr(L)
dat = []
for line in open('data/poly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  dat.append(float(temp[1]))

# Print number of measurements to allow offline checks
print "%d measurements for %s" % (len(dat), tag)

# Create histogram with 20 bins set by the arange
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, bins=np.arange(0.0, 1.0, 0.05),
         log=False, normed=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, 1.0])
plt.grid(False)

plt.xlabel('|PL|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
outfile = 'PL_hist_' + tag + '.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
