#!/usr/bin/python
import sys
import time
from numpy import sqrt, sin, pi, dot, exp
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement four-vector
# Use finite-volume discrete Fourier transform

# Parse arguments: Four-component integer vector and L
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <n_t> <L>"
  sys.exit(1)
n_x = int(sys.argv[1])
n_y = int(sys.argv[2])
n_z = int(sys.argv[3])
n_t = int(sys.argv[4])
L = int(sys.argv[5])
runtime = -time.time()

# Summation range will be -L / 2 + 1, ..., L / 2 inclusive
if not L % 2 == 0:
  print "Error: Need even L rather than", L
  sys.exit(1)
low  = 1 - L / 2
high = L / 2 + 1    # Account for range() dropping upper limit

# For later convenience
invSqrt2  = 2.0 * pi / float(L * sqrt(2.0))
invSqrt6  = 2.0 * pi / float(L * sqrt(6.0))
invSqrt12 = 2.0 * pi / float(L * sqrt(12.0))
halfSqrt5 = pi * sqrt(5.0) / float(L)

# Convert n_i to r_i using usual ehat basis
tag = [n_x - n_y,                 n_x + n_y - 2 * n_z,
       n_x + n_y + n_z - 3 * n_t, n_x + n_y + n_z + n_t]
print("tag: %d, %d, %d, %d" % (tag[0], tag[1], tag[2], tag[3]))
r = [tag[0] / sqrt(2.0),  tag[1] / sqrt(6.0),
     tag[2] / sqrt(12.0), tag[3] / sqrt(20.0)]
mag = sqrt(dot(r, r))
print "|r| = %.4g" % mag
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (p1, p2, p3), each from 0 to L-1,
# except for (0, 0, 0), which must be treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
toSum = 0.0 + 0.0j
for p1 in range(low, high):
  for p2 in range(low, high):
    for p3 in range(low, high):
      for p4 in range(low, high):
        # Omit zero-mode contribution
        # TODO: Separate analytical computation possible?
        if p1 == 0 and p2 == 0 and p3 == 0 and p4 == 0:
          continue

        # Convert p_i to k_i using ghat basis
        # Pattern resembles tag above, now with 2pi factors
        k = [0.0, 0.0, 0.0, 0.0]
        k[0] = (p1 - p2) * invSqrt2
        k[1] = (p1 + p2 - 2.0 * p3) * invSqrt6
        k[2] = (p1 + p2 + p3 - 3.0 * p4) * invSqrt12
        k[3] = (p1 + p2 + p3 + p4) * halfSqrt5

        # Accumulate exp(i r.k) / khatSq (overall factor of 4 outside loop)
        num = exp(1.0j * (r[0] * k[0] + r[1] * k[1] + r[2] * k[2] \
                                                    + r[3] * k[3]))
        denom = (sin(0.5 * k[0]))**2 + (sin(0.5 * k[1]))**2 \
              + (sin(0.5 * k[2]))**2 + (sin(0.5 * k[3]))**2
        toSum += num / denom
#        print "%d %d %d --> %.4g / %.4g = %.4g" \
#              % (p1, p2, p3, num, khatSq, num / khatSq)

# Overall square root and constant factor of sqrt(4pi^2 / 4L^4),
# Determinants now cancel, so include only factor of 2 for trace normalization
one_ov_rI = pi * sqrt(2.0 * toSum) / float(L**2)

# Print along with r_I itself
rI = 1.0 / (one_ov_rI)
print "(%.4g, %.4g) --> (%.4g, %.4g)" \
      % (one_ov_rI.real, one_ov_rI.imag, rI.real, rI.imag)

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
