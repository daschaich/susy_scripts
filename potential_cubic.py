#!/usr/bin/python
import sys
import time
from numpy import sqrt, sin, exp, pi
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

# For later convenience
twopiOvL  = 2.0 * pi / float(L)

# Print naive |r|
mag = sqrt(n_x * n_x + n_y * n_y + n_z * n_z)
print "|r| = %.4g" % mag
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (p1, p2, p3), each from 0 to L-1,
# except for (0, 0, 0), which must be treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
one_ov_rI = 0.0
for p1 in range(L):
  for p2 in range(L):
    for p3 in range(L):
      # Omit zero-mode contribution
      # TODO: Separate analytical computation possible?
      if p1 == 0 and p2 == 0 and p3 == 0:
        continue

      # Accumulate exp(i r.k) / khatSq
      k = [twopiOvL * p1, twopiOvL * p2, twopiOvL * p3]
      num = exp(1.0j * (n_x * k[0] + n_y * k[1] + n_z * k[2]))
      denom = (sin(0.5 * k[0]))**2 + (sin(0.5 * k[1]))**2 \
                                   + (sin(0.5 * k[2]))**2
      one_ov_rI += num / denom
#      print "%d %d %d --> %.4g / %.4g = %.4g" \
#            % (p1, p2, p3, num, khatSq, num / khatSq)

# Constant overall factor of 4pi / 4L^3
# Use standard Tr[T^A T^B] = 0.5\de^{AB} trace normalization...
one_ov_rI *= pi / float(L**3)

# Print along with r_I itself
rI = 1.0 / (one_ov_rI)
print "(%.4g, %.4g) --> (%.4g, %.4g)" \
      % (one_ov_rI.real, one_ov_rI.imag, rI.real, rI.imag)

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
