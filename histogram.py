#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
# ------------------------------------------------------------------
# Plot histogram of data in three files:
#   plaq_$tag.dat, det_re_$tag.dat and det_im_$tag.dat
# Construct these files to include only thermalized data
# Save resulting plot as hist_$tag.pdf
# Adapted from topological charge script from Meifeng

# Parse argument: the tag specifying which files to analyze
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<tag>"
  sys.exit(1)
tag = str(sys.argv[1])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read the data, finding the extremal values
for i in ['plaq_', 'det_re_', 'det_im_']:
  filename = i + tag + '.dat'
  if not os.path.isfile(filename):
    print "ERROR:", filename, "does not exist"
    sys.exit(1)

  if 'plaq' in i:
    plaq_dat = np.loadtxt(filename, unpack=True)
  elif 're' in i:
    re_dat = np.loadtxt(filename, unpack=True)
  elif 'im' in i:
    im_dat = np.loadtxt(filename, unpack=True)

# Determine appropriate order of magnitude for horizontal axis
if max(plaq_dat) - 25.0 > min(plaq_dat):
  base = int(10)
else:
  base = int(1)
maxPlaq = base * int(np.ceil(max(plaq_dat) / float(base)))
minPlaq = base * int(np.floor(min(plaq_dat) / float(base)))
maxRe = base * int(np.ceil(max(re_dat) / float(base)))
minRe = base * int(np.floor(min(re_dat) / float(base)))
maxIm = base * int(np.ceil(max(im_dat) / float(base)))
minIm = base * int(np.floor(min(im_dat) / float(base)))

# Figure out range necessary to combine real and imaginary parts of det
if maxIm > maxRe: maxDet = maxIm
else:             maxDet = maxRe
if minIm < minRe: minDet = minIm
else:             minDet = minRe

# Create an instance of plt and position the subplots
# (MNi) denotes the ith plot in an MxN array
fig = plt.figure()
histogram_plaq = fig.add_subplot(211, xlim = (minPlaq, maxPlaq))
histogram_det = fig.add_subplot(212, xlim = (minDet, maxDet))
fig.subplots_adjust(hspace = 0.35)

# Plaquette histogram plot
nbins = 30
skip = base * int(round((maxPlaq - minPlaq) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_plaq.hist(plaq_dat[:], nbins, log=True, normed=True, align='mid',
                    facecolor='blue', alpha=0.9)

histogram_plaq.set_title(tag)
histogram_plaq.set_xlabel('plaq')
histogram_plaq.set_ylabel('Frequency')
histogram_plaq.set_xticks(range(minPlaq, maxPlaq + 1, skip))

# Plaquette determinant histogram plot
nbins = int(round(30 * (maxRe - minRe) / (maxDet - minDet)))
histogram_det.hist(re_dat, nbins, log=True, normed=True, align='mid',
                   facecolor='green', alpha=0.5, label='Real')

nbins = int(round(30 * (maxIm - minIm) / (maxDet - minDet)))
histogram_det.hist(im_dat, nbins, log=True, normed=True, align='mid',
                   facecolor='red', alpha=0.5, label='Imag')

histogram_det.set_xlabel('det')
histogram_det.set_ylabel('Frequency')
skip = base * int(round((maxDet - minDet) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_det.set_xticks(range(minDet, maxDet + 1, skip))
histogram_det.legend()

# Save a pdf
outfile = 'distribution_' + tag + '.pdf'
plt.savefig(outfile)
# ------------------------------------------------------------------
