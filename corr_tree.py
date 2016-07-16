#!/usr/bin/python
import sys
import time
from numpy import sqrt, sin, exp, pi, dot
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement four-vector
# Use finite-volume discrete Fourier transform

# Parse arguments: Four-component integer vector and L
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<n_1> <n_2> <n_3> <n_4> <L>"
  sys.exit(1)
n = [0, 0, 0, 0]
n[0] = int(sys.argv[1])
n[1] = int(sys.argv[2])
n[2] = int(sys.argv[3])
n[3] = int(sys.argv[4])
L = int(sys.argv[5])
runtime = -time.time()

# A4* basis vectors as 4x4 matrix
invSqrt2  = 1.0 / sqrt(2.0)
invSqrt6  = 1.0 / sqrt(6.0)
invSqrt12 = 1.0 / sqrt(12.0)
invSqrt20 = 1.0 / sqrt(20.0)
pi_ov_L   = pi / float(L)
twopiOvL  = 2.0j * pi / float(L)
ehat = [[invSqrt2,        invSqrt6,        invSqrt12,        invSqrt20], \
        [-1.0 * invSqrt2, invSqrt6,        invSqrt12,        invSqrt20], \
        [0.0,             -2.0 * invSqrt6, invSqrt12,        invSqrt20], \
        [0.0,             0.0,             -3.0 * invSqrt12, invSqrt20]]

# Print naive |r| using ehat
r = [0.0, 0.0, 0.0, 0.0]
for i in range(4):
  for j in range(4):
    r[j] += n[i] * ehat[i][j]

mag = sqrt(dot(r, r))
print "|r| = %.4g" % mag
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (p1, p2, p3, p4), each from 0 to L-1,
# except for (0, 0, 0, 0), which must be treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
toSum = 0.0 + 0.0j
for p1 in range(L):
  for p2 in range(L):
    for p3 in range(L):
      for p4 in range(L):
        # Omit zero-mode contribution
        # TODO: Separate analytical computation possible?
        if p1 + p2 + p3 + p4 == 0:
          continue

        # Accumulate exp(i r.k) / khatSq (overall factor of 4 outside loop)
        num = exp(twopiOvL * (n[0] * p1 + n[1] * p2 + n[2] * p3 + n[3] * p4))
        denom = (sin(pi_ov_L * p1))**2 + (sin(pi_ov_L * p2))**2 \
              + (sin(pi_ov_L * p3))**2 + (sin(pi_ov_L * p4))**2
        toSum += num / denom

# Overall square root and constant factor of sqrt(4pi^2 / 4L^4),
# Don't add any extra normalization (determinants cancel)
one_ov_rI = pi * sqrt(toSum) / float(L**2)

# Print along with r_I itself
rI = 1.0 / (one_ov_rI)
print "(%.4g, %.4g) --> (%.4g, %.4g)" \
      % (one_ov_rI.real, one_ov_rI.imag, rI.real, rI.imag)

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
