#!/usr/bin/python
import sys
import time
import numpy as np
# ------------------------------------------------------------------
# Print r_I for simple cubic lattice as check
# Use finite-volume discrete Fourier transform

# Parse arguments: Three-component integer vector and L
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <L>"
  sys.exit(1)
n_x = int(sys.argv[1])
n_y = int(sys.argv[2])
n_z = int(sys.argv[3])
L = int(sys.argv[4])
runtime = -time.time()

# Summation range will be -L / 2 + 1, ..., L / 2 inclusive
if not L % 2 == 0:
  print "Error: Need even L rather than", L
  sys.exit(1)
low  = 1 - L / 2
high = L / 2 + 1    # Account for range() dropping upper limit

# For later convenience
twopiOvL  = 2.0 * np.pi / float(L)

# Print naive |r|
mag = np.sqrt(n_x * n_x + n_y * n_y + n_z * n_z)
print "|r| = %.4g" % mag
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (p1, p2, p3), each from 0 to L-1,
# except for (0, 0, 0), which is treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
one_ov_rI = 0.0
for p1 in range(low, high):
  for p2 in range(low, high):
    for p3 in range(low, high):
      # Omit zero-mode contribution
      # TODO: Separate analytical computation possible?
      if p1 == 0 and p2 == 0 and p3 == 0:
        continue

      # Convert p_i to k_i using ghat basis
      # Pattern is same as tag above, but now have overall twopiOvL factor
      temp = np.array([p1, p2, p3], dtype = np.float)
      k = twopiOvL * temp
      khat_mu = 2.0 * np.sin(0.5 * k)
      khatSq = (khat_mu**2).sum()

      # Sum \prod_i cos(r_i k_i) / khatSq
      num = np.cos(n_x * k[0]) * np.cos(n_y * k[1]) * np.cos(n_z * k[2])
      one_ov_rI += num / khatSq
#      print "%d %d %d --> %.4g / %.4g = %.4g" \
#            % (p1, p2, p3, num, khatSq, num / khatSq)

# Constant overall factor of 4pi / L^3
# Use standard Tr[T^A T^B] = 0.5\de^{AB} trace normalization...
one_ov_rI *= 4.0 * np.pi / float(L**3)

# Print along with r_I itself
rI = 1.0 / (one_ov_rI)
print "%.4g --> %.4g" % (one_ov_rI, rI)

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
