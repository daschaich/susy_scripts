#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compare 6nt6 and 12nt12 Wilson loops,
# with and without modifications for R symmetry tests

# No arguments yet...
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
for i in ['0.5', '1.0', '1.5']:
  smallFile = 'Nc2_6nt6/l' + i + '_b0.5_f0.0_k0.5/results/rsymm.dat'
  largeFile = 'Nc2_12nt12/l' + i + '_b0.5_f0.0_k0.5/results/rsymm.dat'
  if not os.path.isfile(smallFile) or not os.path.isfile(largeFile):
    print "ERROR: missing file for lambda =", i
    sys.exit(1)

  # Set rescaling parameter xi from 1x1 vs. 2x2 Wilson loops
  # Regular Wilson loops are last two columns
  # Modified Wilson loops are the two before that
  for line in open(smallFile):
    if line.startswith('1 1 '):
      temp = line.split()
      smallPlaq = float(temp[8])
      smallErr = float(temp[9]) / smallPlaq
      break                       # Done with file for now
  for line in open(largeFile):
    if line.startswith('2 2 '):
      temp = line.split()
      largePlaq = float(temp[8])
      largeErr = float(temp[9]) / largePlaq
      break                       # Done with file for now
  xi = smallPlaq / largePlaq      # Technically xi^4
  xi_err = xi * np.sqrt(smallErr**2 + largeErr**2)
#  print "1 1 %.4g %.6g %.4g" % (float(i), xi, xi_err)

  # Check other Wilson loops on 6nt6 and 12nt12
  for line in open(smallFile):
    if line.startswith('1 1 ') or line.startswith('# '):
      continue

    temp = line.split()
    norm = 2 * int(temp[0])
    inv = 2 * int(temp[1])
    if inv < norm:   # Ordinary Wilson loops are symmetric
      continue
    tag = str(norm) + ' ' + str(inv) + ' '
    smallLoop = float(temp[8])
    smallErr = float(temp[9])

    for inner_line in open(largeFile):
      if inner_line.startswith(tag):
        temp = inner_line.split()
        largeLoop = float(temp[8])
        largeErr = float(temp[9])
        break                       # Done with inner loop

    p = float(norm + inv) / 4.0
    rescaled = largeLoop * np.power(xi, p)   # de(xi^p) = p (xi^{p - 1}) de(x)
    err1 = largeErr / largeLoop
    err2 = (p * np.power(xi, p - 1.0) * xi_err) / xi
    err = rescaled * np.sqrt(err1**2 + err2**2)
#    print "%d %d %.4g" % (int(norm) / 2, int(inv) / 2, float(i)),
#    print "%.6g %.4g %.6g %.4g %.6g %.4g" \
#          % (smallLoop, smallErr, largeLoop, largeErr, rescaled, err)

  # Now check blocked rescaled modified Wilson loops
  for line in open(smallFile):
    if line.startswith('# '):
      continue

    temp = line.split()
    norm = 2 * int(temp[0])
    inv = 2 * int(temp[1])
    tag = str(norm) + ' ' + str(inv) + ' '
    smallRel = float(temp[4])
    smallErr = float(temp[5])
    for inner_line in open(largeFile):
      if inner_line.startswith(tag):
        temp = inner_line.split()
        largeRel = float(temp[4])
        largeErr = float(temp[5])
        Mloop = float(temp[6])
        Mtemp = float(temp[7]) / Mloop
        Wloop = float(temp[8])
        Wtemp = float(temp[9]) / Wloop
        break                       # Done with inner loop

    # Rescale ordinary Wilson loop
    p = float(norm + inv) / 4.0
    Wloop *= np.power(xi, p)
    xi_err = (p * np.power(xi, p - 1.0) * xi_err) / xi
    Werr = rescaled * np.sqrt(xi_err**2 + Wtemp**2)

    # Rescale modified Wilson loop
    p = float(norm - inv) / 4.0
    Mloop *= np.power(xi, p)
    xi_err = (p * np.power(xi, p - 1.0) * xi_err) / xi
    Merr = rescaled * np.sqrt(xi_err**2 + Mtemp**2)

    # Compute rescaled (M - W) / N(M + W), propagating uncertainties
    num = (Mloop - Wloop) / inv
    num_err = np.sqrt(Merr**2 + Werr**2) / (inv * num)
    den = Mloop + Wloop
    den_err = np.sqrt(Merr**2 + Werr**2) / den
    ratio = num / den
    err = ratio * np.sqrt(num_err**2 + den_err**2)

    # Print out in same format as rsymm.dat
#    print "%d %d" % (int(norm) / 2, int(inv) / 2),
#    print "%.6g %.4g" % (num, num_err),
#    print "%.6g %.4g" % (ratio, err),
#    print "%.6g %.4g" % (Mloop, Merr),
#    print "%.6g %.4g" % (Wloop, Werr)

    print "%d %d %.4g" % (int(norm) / 2, int(inv) / 2, float(i)),
    print "%.6g %.4g %.6g %.4g %.6g %.4g" \
          % (smallRel, smallErr, largeRel, largeErr, ratio, err)
# ------------------------------------------------------------------
