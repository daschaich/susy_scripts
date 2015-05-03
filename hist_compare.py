#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of data in ten files:
#   plaq_$tag.dat, tr_$tag.dat, det_$tag.dat,
#   re_$tag.dat and im_$tag.dat
# for each of two input tags
# Construct these files to include only thermalized data
# Save resulting plot as comparison$tag.pdf
# Adapted from topological charge script from Meifeng

# Parse argument: the tags specifying which files to analyze
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<tag1> <tag2>"
  sys.exit(1)
tag = str(sys.argv[1])
comp = str(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read the data, finding the extremal values
for i in ['plaq_', 'tr_', 'det_', 're_', 'im_']:
  filename = i + tag + '.dat'
  if not os.path.isfile(filename):
    print "ERROR:", filename, "does not exist"
    sys.exit(1)

  if 'plaq' in i:
    plaq_dat = np.loadtxt(filename, unpack=True)
  elif 'tr' in i:
    tr_dat = np.loadtxt(filename, unpack=True)
  elif 'det' in i:
    det_dat = np.loadtxt(filename, unpack=True)
  elif 're' in i:
    re_dat = np.loadtxt(filename, unpack=True)
  elif 'im' in i:
    im_dat = np.loadtxt(filename, unpack=True)

  # Other file for comparison
  filename = i + comp + '.dat'
  if not os.path.isfile(filename):
    print "ERROR:", filename, "does not exist"
    sys.exit(1)

  if 'plaq' in i:
    plaq_comp = np.loadtxt(filename, unpack=True)
  elif 'tr' in i:
    tr_comp = np.loadtxt(filename, unpack=True)
  elif 'det' in i:
    det_comp = np.loadtxt(filename, unpack=True)
  elif 're' in i:
    re_comp = np.loadtxt(filename, unpack=True)
  elif 'im' in i:
    im_comp = np.loadtxt(filename, unpack=True)

# Determine appropriate order of magnitude for horizontal axes
plaq_all = np.concatenate((plaq_dat, plaq_comp), axis=0)
tr_all = np.concatenate((tr_dat, tr_comp), axis=0)
det_all = np.concatenate((det_dat, det_comp), axis=0)
re_all = np.concatenate((re_dat, re_comp), axis=0)
im_all = np.concatenate((im_dat, im_comp), axis=0)
if max(plaq_all) - 25.0 > min(plaq_all):
  base = int(10)
else:
  base = int(1)
maxPlaq = base * int(np.ceil(max(plaq_all) / float(base)))
minPlaq = base * int(np.floor(min(plaq_all) / float(base)))
maxTr = base * int(np.ceil(max(tr_all) / float(base)))
minTr = base * int(np.floor(min(tr_all) / float(base)))
maxDet = base * int(np.ceil(max(det_all) / float(base)))
minDet = base * int(np.floor(min(det_all) / float(base)))
maxRe = base * int(np.ceil(max(re_all) / float(base)))
minRe = base * int(np.floor(min(re_all) / float(base)))
maxIm = base * int(np.ceil(max(im_all) / float(base)))
minIm = base * int(np.floor(min(im_all) / float(base)))

# Create an instance of plt and position the subplots
# (MNi) denotes the ith plot in an MxN array
# Stretch from default 8"x6" to 8"x15" for five figures
fig = plt.figure(figsize = (8, 15))
histogram_plaq = fig.add_subplot(511, xlim = (minPlaq, maxPlaq))
histogram_tr = fig.add_subplot(512, xlim = (minTr, maxTr))
histogram_det = fig.add_subplot(513, xlim = (minDet, maxDet))
histogram_re = fig.add_subplot(514, xlim = (minRe, maxRe))
histogram_im = fig.add_subplot(515, xlim = (minIm, maxIm))
fig.subplots_adjust(hspace = 0.35)

# Plaquette histogram plot
nbins = 30
histogram_plaq.hist(plaq_dat, nbins, log=True, normed=True, align='mid',
                    facecolor='blue', alpha=0.5, label=tag,
                    histtype='stepfilled')
histogram_plaq.hist(plaq_comp, nbins, log=True, normed=True, align='mid',
                    facecolor='red', alpha=0.5, label=comp,
                    histtype='stepfilled')
title = ((os.getcwd()).split('/'))[-1]
histogram_plaq.set_title(title)
histogram_plaq.set_xlabel('plaq')
histogram_plaq.set_ylabel('Frequency')
skip = base * int(round((maxPlaq - minPlaq) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_plaq.set_xticks(range(minPlaq, maxPlaq + 1, skip))
histogram_plaq.legend()

# Scalar potential histogram plot
histogram_tr.hist(tr_dat, nbins, log=True, normed=True, align='mid',
                  facecolor='blue', alpha=0.5, label=tag,
                  histtype='stepfilled')
histogram_tr.hist(tr_comp, nbins, log=True, normed=True, align='mid',
                  facecolor='red', alpha=0.5, label=comp,
                  histtype='stepfilled')
histogram_tr.set_xlabel('Tr[U.Udag]/N-1')
histogram_tr.set_ylabel('Frequency')
skip = base * int(round((maxTr - minTr) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_tr.set_xticks(range(minTr, maxTr + 1, skip))
histogram_tr.legend()

# Plaquette determinant histogram plot
histogram_det.hist(det_dat, nbins, log=True, normed=True, align='mid',
                  facecolor='blue', alpha=0.5, label=tag,
                  histtype='stepfilled')
histogram_det.hist(det_comp, nbins, log=True, normed=True, align='mid',
                  facecolor='red', alpha=0.5, label=comp,
                  histtype='stepfilled')
histogram_det.set_xlabel('|det - 1|^2')
histogram_det.set_ylabel('Frequency')
skip = base * int(round((maxDet - minDet) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_det.set_xticks(range(minDet, maxDet + 1, skip))
histogram_det.legend()

# Plaquette determinant real and imaginary parts
histogram_re.hist(re_dat, nbins, log=True, normed=True, align='mid',
                  facecolor='blue', alpha=0.5, label=tag,
                  histtype='stepfilled')
histogram_re.hist(re_comp, nbins, log=True, normed=True, align='mid',
                  facecolor='red', alpha=0.5, label=comp,
                  histtype='stepfilled')
histogram_re.set_xlabel('Re(det)')
histogram_re.set_ylabel('Frequency')
skip = base * int(round((maxRe - minRe) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_re.set_xticks(range(minRe, maxRe + 1, skip))
histogram_re.legend()

histogram_im.hist(im_dat, nbins, log=True, normed=True, align='mid',
                  facecolor='blue', alpha=0.5, label=tag,
                  histtype='stepfilled')
histogram_im.hist(im_comp, nbins, log=True, normed=True, align='mid',
                  facecolor='red', alpha=0.5, label=comp,
                  histtype='stepfilled')

histogram_im.set_xlabel('Im(det)')
histogram_im.set_ylabel('Frequency')
skip = base * int(round((maxIm - minIm) / (8.0 * float(base))))
if skip == 0:
  skip = 1
histogram_im.set_xticks(range(minIm, maxIm + 1, skip))
histogram_im.legend()

# Save a pdf
outfile = 'comparison_' + tag + '.pdf'
plt.savefig(outfile)
# ------------------------------------------------------------------
