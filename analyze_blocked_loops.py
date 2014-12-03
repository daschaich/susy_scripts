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
#blmax = 2
#for la in ['0.5', '1.0', '1.5']:
#  smallFile = 'Nc2_6nt6/l' + i + '_b0.5_f0.0_k0.5/results/rsymm.dat'
#  largeFile = 'Nc2_12nt12/l' + i + '_b0.5_f0.0_k0.5/results/rsymm.dat'
blmax = 3
for la in ['1.0']:
  smallFile = 'Nc2_8nt16/l' + la + '_b1.0_f0.0_k1.0/results/rsymm.dat'
  largeFile = 'Nc2_16nt32/l' + la + '_b1.0_f0.0_k1.0/results/rsymm.dat'
  if not os.path.isfile(smallFile) or not os.path.isfile(largeFile):
    print "ERROR: missing file for lambda =", la
    sys.exit(1)

  # Print out rescaling parameter xi for each blocking level
  # Require that n-times-blocked (2n)x(2n) Wilson loop
  # matches (n-1)-times-blocked nxn Wilson loop
  # Regular Wilson loops are last two columns
  # Modified Wilson loops are the two before that
  xi = np.zeros(blmax, dtype = np.float)
  xiErr = np.zeros(blmax, dtype = np.float)
  smallPlaq = np.zeros(blmax, dtype = np.float)
  smallErr = np.zeros(blmax, dtype = np.float)
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    norm = int(temp[0])
    inv  = int(temp[1])
    if not norm == inv:
      continue
    if norm == 1:
      smallPlaq[0] = float(temp[8])
      smallErr[0] = float(temp[9]) / float(temp[8])
    elif norm == 2:
      smallPlaq[1] = float(temp[8])
      smallErr[1] = float(temp[9]) / float(temp[8])
    elif norm == 4:
      smallPlaq[2] = float(temp[8])
      smallErr[2] = float(temp[9]) / float(temp[8])

  largePlaq = np.zeros(blmax, dtype = np.float)
  largeErr = np.zeros(blmax, dtype = np.float)
  for bl in range(blmax - 1):
    largePlaq[bl] = smallPlaq[bl + 1]
    largeErr[bl] = smallErr[bl + 1]
  for line in open(largeFile):
    if line.startswith('8 8 '):
      temp = line.split()
      largePlaq[-1] = float(temp[8])
      largeErr[-1] = float(temp[9]) / float(temp[8])

  print "0 1.0     0.0"
  for bl in range(blmax):
    xi[bl] = smallPlaq[bl] / largePlaq[bl]      # Technically xi^4
    xiErr[bl] = xi[bl] * np.sqrt(smallErr[bl]**2 + largeErr[bl]**2)
    print "%d %.6g %.4g" % (bl + 1, xi[bl], xiErr[bl])
  sys.exit(0)

  # Check other Wilson loops
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
