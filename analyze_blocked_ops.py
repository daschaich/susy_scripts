#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Determine xi from small-volume and blocked-large-volume Konishi operators,
# scaling mu \propto 1 / L

# Only consider operator built from two log-polar links
# Equating these should give xi^2; print out xi

# No arguments yet...
small_tag = []
large_tag = []
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# U(4) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 4
#large_dir = 'Nc4_12nt12/'
#small_dir = 'Nc4_6nt6/'
#large_tag.append('l1.0_b0.15_G0.05');  small_tag.append('l1.0_b0.3_G0.05');

# U(4) 8nt8 & 4nt4
#Nc = 4
#large_dir = 'Nc4_8nt8/'
#small_dir = 'Nc4_4nt4/'
#large_tag.append('l1.0_b0.2_G0.05');    small_tag.append('l1.0_b0.4_G0.05');
#large_tag.append('l2.0_b0.3_G0.05');    small_tag.append('l2.0_b0.57_G0.05');

# U(3) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 3
#large_dir = 'Nc3_12nt12/'
#small_dir = 'Nc3_6nt6/'
#large_tag.append('l1.0_b0.15_G0.05');   small_tag.append('l1.0_b0.27_G0.05')
#large_tag.append('l2.0_b0.2_G0.05');    small_tag.append('l2.0_b0.4_G0.05')
#large_tag.append('l3.0_b0.23_G0.05');   small_tag.append('l3.0_b0.5_G0.05')

# U(3) 8nt8 & 4nt4 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 3
#large_dir = 'Nc3_8nt8/'
#small_dir = 'Nc3_4nt4/'
#large_tag.append('l1.0_b0.2_G0.05');    small_tag.append('l1.0_b0.4_G0.05')
#large_tag.append('l2.0_b0.3_G0.05');    small_tag.append('l2.0_b0.57_G0.05')
#large_tag.append('l3.0_b0.3_G0.05');    small_tag.append('l3.0_b0.69_G0.05') # Runs away...

# U(2) 16nt16 & 8nt8 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 2
#large_dir = 'Nc2_16nt16/'
#small_dir = 'Nc2_8nt8/'
#large_tag.append('l0.5_b0.07_G0.05');   small_tag.append('l0.5_b0.14_G0.05')
#large_tag.append('l1.0_b0.1_G0.05');    small_tag.append('l1.0_b0.2_G0.05')
#large_tag.append('l2.0_b0.14_G0.05');   small_tag.append('l2.0_b0.28_G0.05')
#large_tag.append('l3.0_b0.17_G0.05');   small_tag.append('l3.0_b0.35_G0.05')
#large_tag.append('l4.0_b0.2_G0.05');    small_tag.append('l4.0_b0.4_G0.05')
#large_tag.append('l5.0_b0.22_G0.05');   small_tag.append('l5.0_b0.44_G0.05')

# U(2) 12nt12 & 6nt6 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 2
#large_dir = 'Nc2_12nt12/'
#small_dir = 'Nc2_6nt6/'
#large_tag.append('l0.5_b0.095_G0.05');  small_tag.append('l0.5_b0.19_G0.05')
#large_tag.append('l1.0_b0.13_G0.05');   small_tag.append('l1.0_b0.26_G0.05')
#large_tag.append('l2.0_b0.19_G0.05');   small_tag.append('l2.0_b0.38_G0.05')
#large_tag.append('l3.0_b0.23_G0.05');   small_tag.append('l3.0_b0.46_G0.05')
#large_tag.append('l4.0_b0.25_G0.05');   small_tag.append('l4.0_b0.5_G0.05')
#large_tag.append('l5.0_b0.3_G0.05');    small_tag.append('l5.0_b0.6_G0.05')

# U(2) 8nt8 & 4nt4 with G=0.05 and (mu*L)^2 / lambda ~ 2.5
#Nc = 2
#large_dir = 'Nc2_8nt8/'
#small_dir = 'Nc2_4nt4/'
#large_tag.append('l0.5_b0.14_G0.05');   small_tag.append('l0.5_b0.28_G0.05')
#large_tag.append('l1.0_b0.2_G0.05');    small_tag.append('l1.0_b0.4_G0.05')
#large_tag.append('l2.0_b0.28_G0.05');   small_tag.append('l2.0_b0.57_G0.05')
#large_tag.append('l3.0_b0.35_G0.05');   small_tag.append('l3.0_b0.69_G0.05')
#large_tag.append('l4.0_b0.4_G0.05');    small_tag.append('l4.0_b0.8_G0.05')

for i in range(len(small_tag)):
  outfilename = large_dir + large_tag[i] + '/results/xi_op.dat'
  print "Comparing %s vs. %s" % (small_dir + small_tag[i], \
                                 large_dir + large_tag[i])
  smallFile = small_dir + small_tag[i] + '/results/ops.dat'
  largeFile = large_dir + large_tag[i] + '/results/ops.dat'

  # First make sure we're calling this from the right place
  # Should be able to retain these independent of commenting out above
  if not os.path.isfile(smallFile):
    print "ERROR: missing file", smallFile
    sys.exit(1)
  if not os.path.isfile(largeFile):
    print "ERROR: missing file", largeFile
    sys.exit(1)

  # Print out rescaling parameter xi for each blocking level
  # Require that n-times-blocked L Konishi operator
  # matches (n-1)-times-blocked (L/2) result
  smallOp  = []
  smallErr = []
  for line in open(smallFile):
    if line.startswith('# '):
      continue
    temp = line.split()
    smallOp.append(float(temp[1]))
    smallErr.append(float(temp[2]) / float(temp[1]))
  largeOp  = []
  largeErr = []
  for line in open(largeFile):
    if line.startswith('# ') or line.startswith('0 '):
      continue
    temp = line.split()
    largeOp.append(float(temp[1]))
    largeErr.append(float(temp[2]) / float(temp[1]))

  blmax = len(smallOp)
  if not len(largeOp) == blmax:
    print "ERROR: inconsistent lengths: %d %d" % (blmax, len(largeOp))
    sys.exit(1)

  outfile = open(outfilename, 'w')
  print "0 1.0     0.0"
  print >> outfile, "0 1.0     0.0"
  for bl in range(blmax):
    xi = np.sqrt(smallOp[bl] / largeOp[bl])   # Square root to get xi itself
                                              # -->0.5 in error propagation
    xi_err = 0.5 * xi * np.sqrt(smallErr[bl]**2 + largeErr[bl]**2)
    print "%d %.6g %.4g" % (bl + 1, xi, xi_err)
    print >> outfile, "%d %.6g %.4g" % (bl + 1, xi, xi_err)
  outfile.close()
# ------------------------------------------------------------------
