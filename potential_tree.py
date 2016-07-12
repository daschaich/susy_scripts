#!/usr/bin/python
import sys
import time
from numpy import sqrt, sin, exp, pi, dot
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement three-vector
# Use finite-volume discrete Fourier transform

# Parse arguments: Three-component integer vector and L
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<n_1> <n_2> <n_3> <L>"
  sys.exit(1)
n = [0, 0, 0]
n[0] = int(sys.argv[1])
n[1] = int(sys.argv[2])
n[2] = int(sys.argv[3])
L = int(sys.argv[4])
runtime = -time.time()

# Timeslice basis vectors as 3x3 matrix
invSqrt2  = 1.0 / sqrt(2.0)
invSqrt6  = 1.0 / sqrt(6.0)
invSqrt12 = 1.0 / sqrt(12.0)
pi_ov_L   = pi / float(L)
twopiOvL  = 2.0j * pi / float(L)
ehat = [[invSqrt2,        invSqrt6,        invSqrt12], \
        [-1.0 * invSqrt2, invSqrt6,        invSqrt12], \
        [0.0,             -2.0 * invSqrt6, invSqrt12]]

# Print naive |r| using ehat
r = [0.0, 0.0, 0.0]
for i in range(3):
  for j in range(3):
    r[j] += n[i] * ehat[i][j]

mag = sqrt(dot(r, r))
print "|r| = %.4g" % mag
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (p1, p2, p3), each from 0 to L-1,
# except for (0, 0, 0), which must be treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
one_ov_rI = 0.0 + 0.0j
for p1 in range(L):
  for p2 in range(L):
    for p3 in range(L):
      # Omit zero-mode contribution
      # TODO: Separate analytical computation possible?
      if p1 + p2 + p3 == 0:
        continue

      # Accumulate exp(i r.k) / khatSq
      num = exp(twopiOvL * (n[0] * p1 + n[1] * p2 + n[2] * p3))
      denom = (sin(pi_ov_L * p1))**2 + (sin(pi_ov_L * p2))**2 \
                                     + (sin(pi_ov_L * p3))**2
      one_ov_rI += num / denom

# Constant overall factor of 4pi / 4L^3
# Determinants cancel, so include only factor of 2 for trace normalization
one_ov_rI *= 2.0 * pi / float(L**3)

# Print along with r_I itself
rI = 1.0 / (one_ov_rI)
print "(%.4g, %.4g) --> (%.4g, %.4g)" \
      % (one_ov_rI.real, one_ov_rI.imag, rI.real, rI.imag)

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
