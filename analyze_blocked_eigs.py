#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Determine xi from small-volume and blocked-large-volume scalar eigenvalues,
# scaling mu \propto 1 / L

# No arguments yet...
small_tag = []
large_tag = []
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# U(4) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc4_12nt12/l1.0_b0.15_G0.05');  small_tag.append('Nc4_6nt6/l1.0_b0.3_G0.05')

# U(4) 8nt8 & 4nt4 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc4_8nt8/l1.0_b0.2_G0.05');     small_tag.append('Nc4_4nt4/l1.0_b0.4_G0.05')
large_tag.append('Nc4_8nt8/l2.0_b0.3_G0.05');     small_tag.append('Nc4_4nt4/l2.0_b0.57_G0.05')

# U(3) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc3_12nt12/l1.0_b0.15_G0.05');  small_tag.append('Nc3_6nt6/l1.0_b0.27_G0.05')
large_tag.append('Nc3_12nt12/l2.0_b0.2_G0.05');   small_tag.append('Nc3_6nt6/l2.0_b0.4_G0.05')
large_tag.append('Nc3_12nt12/l3.0_b0.23_G0.05');  small_tag.append('Nc3_6nt6/l3.0_b0.5_G0.05')

# U(3) 8nt8 & 4nt4 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc3_8nt8/l1.0_b0.2_G0.05');     small_tag.append('Nc3_4nt4/l1.0_b0.4_G0.05')
large_tag.append('Nc3_8nt8/l2.0_b0.3_G0.05');     small_tag.append('Nc3_4nt4/l2.0_b0.57_G0.05')

# U(2) 16nt16 & 8nt8 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc2_16nt16/l0.5_b0.07_G0.05');  small_tag.append('Nc2_8nt8/l0.5_b0.14_G0.05')
large_tag.append('Nc2_16nt16/l1.0_b0.1_G0.05');   small_tag.append('Nc2_8nt8/l1.0_b0.2_G0.05')
large_tag.append('Nc2_16nt16/l2.0_b0.14_G0.05');  small_tag.append('Nc2_8nt8/l2.0_b0.28_G0.05')
large_tag.append('Nc2_16nt16/l3.0_b0.17_G0.05');  small_tag.append('Nc2_8nt8/l3.0_b0.35_G0.05')
large_tag.append('Nc2_16nt16/l4.0_b0.2_G0.05');   small_tag.append('Nc2_8nt8/l4.0_b0.4_G0.05')

# U(2) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc2_12nt12/l0.5_b0.095_G0.05'); small_tag.append('Nc2_6nt6/l0.5_b0.19_G0.05')
large_tag.append('Nc2_12nt12/l1.0_b0.13_G0.05');  small_tag.append('Nc2_6nt6/l1.0_b0.26_G0.05')
large_tag.append('Nc2_12nt12/l2.0_b0.19_G0.05');  small_tag.append('Nc2_6nt6/l2.0_b0.38_G0.05')
large_tag.append('Nc2_12nt12/l3.0_b0.23_G0.05');  small_tag.append('Nc2_6nt6/l3.0_b0.46_G0.05')
large_tag.append('Nc2_12nt12/l4.0_b0.25_G0.05');  small_tag.append('Nc2_6nt6/l4.0_b0.5_G0.05')
large_tag.append('Nc2_12nt12/l5.0_b0.3_G0.05');   small_tag.append('Nc2_6nt6/l5.0_b0.6_G0.05')

# U(2) 8nt8 & 4nt4 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
large_tag.append('Nc2_8nt8/l0.5_b0.14_G0.05');    small_tag.append('Nc2_4nt4/l0.5_b0.28_G0.05')
large_tag.append('Nc2_8nt8/l1.0_b0.2_G0.05');     small_tag.append('Nc2_4nt4/l1.0_b0.4_G0.05')
large_tag.append('Nc2_8nt8/l2.0_b0.28_G0.05');    small_tag.append('Nc2_4nt4/l2.0_b0.57_G0.05')
large_tag.append('Nc2_8nt8/l3.0_b0.35_G0.05');    small_tag.append('Nc2_4nt4/l3.0_b0.69_G0.05')
large_tag.append('Nc2_8nt8/l4.0_b0.4_G0.05');     small_tag.append('Nc2_4nt4/l4.0_b0.8_G0.05')

for i in range(len(small_tag)):
  outfilename = large_tag[i] + '/results/xi.eig'
  print "Comparing %s vs. %s" % (small_tag[i], large_tag[i])
  smallFile = small_tag[i] + '/results/blocked_eigs.dat'
  largeFile = large_tag[i] + '/results/blocked_eigs.dat'

  # First make sure we're calling this from the right place
  if not os.path.isfile(smallFile):
    print "ERROR: missing file", smallFile
    sys.exit(1)
  if not os.path.isfile(largeFile):
    print "ERROR: missing file", largeFile
    sys.exit(1)

  # Extract Nc from path
  if large_tag[i].startswith('Nc2'):
    Nc = 2
  elif large_tag[i].startswith('Nc3'):
    Nc = 3
  elif large_tag[i].startswith('Nc4'):
    Nc = 4
  else:
    print "Couldn't determine Nc from", large_tag[i]

  # Print out rescaling parameter xi for each blocking level
  # Require that rescaled n-times-blocked scalar eigenvalues
  # matches the corresponding (n-1)-times-blocked eigenvalue
  smallEig = [[] ,[]]
  smallErr = [[] ,[]]
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    # Format: num bl ave err
    temp = line.split()
    num = int(temp[0])
    if num == 0:
      smallEig[0].append(float(temp[2]))
      smallErr[0].append(float(temp[3]) / float(temp[2]))
    elif num == Nc - 1:
      smallEig[1].append(float(temp[2]))
      smallErr[1].append(float(temp[3]) / float(temp[2]))

  largeEig = [[] ,[]]
  largeErr = [[] ,[]]
  for line in open(largeFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    num = int(temp[0])
    bl = int(temp[1])
    if bl > 0 and num == 0:
      largeEig[0].append(float(temp[2]))
      largeErr[0].append(float(temp[3]) / float(temp[2]))
    elif bl > 0 and num == Nc - 1:
      largeEig[1].append(float(temp[2]))
      largeErr[1].append(float(temp[3]) / float(temp[2]))

  # Sanity check
  blmax = len(smallEig[0])
  if not len(largeEig[0]) == blmax or not len(smallEig[1]) == blmax \
                                   or not len(largeEig[1]) == blmax:
    print "ERROR: inconsistent lengths: %d %d %d %d" \
          % (blmax, len(largeEig[0]), len(smallEig[1]), len(largeEig[1]))
    sys.exit(1)

  # Print xi
  xi = np.zeros((2, blmax), dtype = np.float)
  xiErr = np.zeros((2, blmax), dtype = np.float)
  outfile = open(outfilename, 'w')
  print "0 1.0     0.0"
  print >> outfile, "0 1.0     0.0"
  for eig in range(2):
    for bl in range(blmax):
      xi[eig][bl] = smallEig[eig][bl] / largeEig[eig][bl]
      xiErr[eig][bl] = xi[eig][bl]
      xiErr[eig][bl] *= np.sqrt(smallErr[eig][bl]**2 + largeErr[eig][bl]**2)
  for bl in range(blmax):
    print "%d %.6g %.4g %.6g %.4g" \
          % (bl + 1, xi[0][bl], xiErr[0][bl], xi[1][bl], xiErr[1][bl])
    print >> outfile, "%d %.6g %.4g %.6g %.4g" \
          % (bl + 1, xi[0][bl], xiErr[0][bl], xi[1][bl], xiErr[1][bl])
  print
  outfile.close()
# ------------------------------------------------------------------
