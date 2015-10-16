#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compare small-volume and blocked-large-volume Wilson loops,
# scaling mu \propto 1 / L
# First determine xi^2 from the (blocked) plaquette
# Then use xi^4 to compute larger loop ratios,
# including modified loops to test discrete R symmetries
# Everything after xi^4 printing is currently commented out

# No arguments yet...
small_tag = []
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# U(4) 8nt8 & 4nt4
#small_dir = 'Nc4_4nt4/'
#small_tag = ['l1.0_b0.4_G0.05', 'l2.0_b0.6_G0.05']
#large_dir = 'Nc4_8nt8/'
#large_tag = ['l1.0_b0.2_G0.05', 'l2.0_b0.3_G0.05']

# U(3) 8nt8 & 4nt4
#small_dir = 'Nc3_4nt4/'
#small_tag = ['l1.0_b0.4_G0.05', 'l2.0_b0.6_G0.05', 'l3.0_b0.8_G0.05']
#large_dir = 'Nc3_8nt8/'
#large_tag = ['l1.0_b0.2_G0.05', 'l2.0_b0.3_G0.05', 'l3.0_b0.4_G0.05']

# U(2) 16nt16 & 8nt8
#small_dir = 'Nc2_8nt8/'
#small_tag = ['l1.0_b0.2_G0.05', 'l1.0_b0.4_G0.05', \
#                                'l2.0_b0.4_G0.05']
#large_dir = 'Nc2_16nt16/'
#large_tag = ['l1.0_b0.1_G0.05', 'l1.0_b0.2_G0.05', \
#                                'l2.0_b0.2_G0.05']

# U(2) 12nt12 & 6nt6
#small_dir = 'Nc2_6nt6/'
#small_tag = ['l1.0_b0.25_G0.05', 'l1.0_b0.5_G0.05', \
#                                 'l2.0_b0.5_G0.05', \
#                                 'l3.0_b0.5_G0.05', \
#                                 'l4.0_b0.5_G0.05', \
#                                 'l5.0_b0.5_G0.05']
#large_dir = 'Nc2_12nt12/'
#large_tag = ['l1.0_b0.13_G0.05', 'l1.0_b0.25_G0.05', \
#                                 'l2.0_b0.25_G0.05', \
#                                 'l3.0_b0.25_G0.05', \
#                                 'l4.0_b0.25_G0.05', \
#                                 'l5.0_b0.25_G0.05']

# U(2) 8nt8 & 4nt4
#small_dir = 'Nc2_4nt4/'
#small_tag = ['l0.5_b0.8_G0.1', 'l1.0_b0.4_G0.05', 'l1.0_b0.8_G0.05', \
#                                                  'l2.0_b0.8_G0.05', \
#                                                  'l3.0_b0.8_G0.05', \
#                                                  'l4.0_b0.8_G0.05', \
#             'l5.0_b0.8_G0.1', \
#             'l6.0_b0.8_G0.1']
#large_dir = 'Nc2_8nt8/'
#large_tag = ['l0.5_b0.4_G0.1', 'l1.0_b0.2_G0.05', 'l1.0_b0.4_G0.05', \
#                                                  'l2.0_b0.4_G0.05', \
#                                                  'l3.0_b0.4_G0.05', \
#                                                  'l4.0_b0.4_G0.05',
#             'l5.0_b0.4_G0.1', \
#             'l6.0_b0.4_G0.1']

for i in range(len(small_tag)):
  if small_tag[i] == large_tag[i]:
    outfilename = large_dir + large_tag[i] + '/results/xi_fixed.dat'
    print "Comparing %s for %s vs. %s" % (small_tag[i], small_dir, large_dir)
  else:
    outfilename = large_dir + large_tag[i] + '/results/xi_scaled.dat'
    print "Comparing %s vs. %s" % (small_dir + small_tag[i], \
                                   large_dir + large_tag[i])
  smallFile = small_dir + small_tag[i] + '/results/rsymm.dat'
  largeFile = large_dir + large_tag[i] + '/results/rsymm.dat'

  # First make sure we're calling this from the right place
  # Should be able to retain these independent of commenting out above
  if not os.path.isfile(smallFile):
    print "ERROR: missing file", smallFile
    sys.exit(1)
  if not os.path.isfile(largeFile):
    print "ERROR: missing file", largeFile
    sys.exit(1)

  # Print out rescaling parameter xi for each blocking level
  # Require that n-times-blocked (2n)x(2n) Wilson loop
  # matches (n-1)-times-blocked nxn Wilson loop
  # Regular Wilson loops are last two columns
  # Modified Wilson loops are the two before that
  smallPlaq = []
  smallErr = []
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    norm = int(temp[0])
    inv  = int(temp[1])
    if not norm == inv:
      continue
    if norm == 1 or norm == 2 or norm == 4:
      smallPlaq.append(float(temp[8]))
      smallErr.append(float(temp[9]) / float(temp[8]))

  largePlaq = []
  largeErr = []
  for line in open(largeFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    norm = int(temp[0])
    inv  = int(temp[1])
    if not norm == inv:
      continue
    if norm == 2 or norm == 4 or norm == 8:
      largePlaq.append(float(temp[8]))
      largeErr.append(float(temp[9]) / float(temp[8]))

  # Sanity check
  blmax = len(smallPlaq)
  if not len(largePlaq) == blmax:
    print "ERROR: %s and %s have different lengths: %d and %d" \
          % (smallFile, largeFile, blmax, len(largePlaq))
    sys.exit(1)

  # Print xi
  xi = np.zeros(blmax, dtype = np.float)
  xiErr = np.zeros(blmax, dtype = np.float)
  outfile = open(outfilename, 'w')
#  print "0 1.0     0.0"
  print >> outfile, "0 1.0     0.0"
  for bl in range(blmax):
    xi[bl] = smallPlaq[bl] / largePlaq[bl]      # Technically xi^4
    xiErr[bl] = xi[bl] * np.sqrt(smallErr[bl]**2 + largeErr[bl]**2)
#    print "%d %.6g %.4g" % (bl + 1, xi[bl], xiErr[bl])
    print >> outfile, "%d %.6g %.4g" % (bl + 1, xi[bl], xiErr[bl])
#  print
  outfile.close()
sys.exit(0)
for la in ['0.5', '1.0', '1.5']:

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
