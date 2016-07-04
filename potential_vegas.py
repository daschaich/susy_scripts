#!/usr/bin/python
from __future__ import print_function     # For vegas
import sys
import time
from numpy import sqrt, sin, cos, pi, dot
import vegas
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement three-vector
# Use vegas to numerically calculate the continuum Fourier transform

# Parse arguments: Three-component integer vector
# and number of evaluations in vegas calculation
if len(sys.argv) < 5:
  print("Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <Neval>")
  sys.exit(1)
n_x = int(sys.argv[1])
n_y = int(sys.argv[2])
n_z = int(sys.argv[3])
Neval = int(sys.argv[4])
Nwarm = int(Neval / 10)
runtime = -time.time()

# For later convenience
invSqrt2  = 1.0 / sqrt(2.0)
invSqrt6  = 1.0 / sqrt(6.0)
invSqrt12 = 1.0 / sqrt(12.0)
twopi = 2.0 * pi

# Convert (n_x, n_y, n_z) to (r_1, r_2, r_3)
tag = [n_x - n_y, n_x + n_y - 2 * n_z, n_x + n_y + n_z]
print("tag: %d, %d, %d" % (tag[0], tag[1], tag[2]))
r = [tag[0] * invSqrt2, tag[1] * invSqrt6, tag[2] * invSqrt12]
mag = sqrt(dot(r, r))
print("|r| = %.4g" % mag)

# Function to integrate
# Convert p_i to k_i using ghat basis
# Pattern is same as tag above, but now have overall twopi factor
# Now integrating over dp = dk / (2pi), so no more 2pi factors in measure
# Include factor of 1/2 for determinant squared times trace normalization
def f(p):
  k = [0.0, 0.0, 0.0]
  k[0] = twopi * (p[0] - p[1]) * invSqrt2
  k[1] = twopi * (p[0] + p[1] - 2.0 * p[2]) * invSqrt6
  k[2] = twopi * (p[0] + p[1] + p[2]) * invSqrt12
  num = cos(dot(r, k))                # Can only handle real part of exp
  denom = (sin(0.5 * k[0]))**2 + (sin(0.5 * k[1]))**2 + (sin(0.5 * k[2]))**2
  return 0.5 * pi * num / denom

integ = vegas.Integrator([[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]])
integ(f, nitn=7, neval=Nwarm)             # Initial adaptation
result = integ(f, nitn=20, neval=Neval)   # Actual estimation

#print(result.summary())
print("result = %s    Q = %.2f" % (result, result.Q))

# Print r_I itself with propagated uncertainty delta(1 / r) = delta(r) / r^2
# Some shenanigans necessary to get numerical results out of vegas class
temp = str(result).split('(')
val = 1.0 / float(temp[0])
err = float(((temp[1]).split(')'))[0])
for i in range(len(temp[0]) - 2):   # Restore leading zeroes to error
  err *= 0.1
print("--> r_I = %.8g %.4g" % (val, err * val * val))

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------
