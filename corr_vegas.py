#!/usr/bin/python
from __future__ import print_function     # For vegas
import sys
import time
from numpy import sqrt, sin, exp, pi, dot
import vegas
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement four-vector
# Use vegas to numerically calculate the infinite-volume Fourier transform

# Parse arguments: Four-component integer vector
# and number of evaluations in vegas calculation
if len(sys.argv) < 6:
  print("Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <n_t> <Neval>")
  sys.exit(1)
n_x = int(sys.argv[1])
n_y = int(sys.argv[2])
n_z = int(sys.argv[3])
n_t = int(sys.argv[4])
Neval = int(sys.argv[5])
Nwarm = int(Neval / 10)
runtime = -time.time()

# For later convenience
invSqrt2  = 2.0 * pi / sqrt(2.0)
halfSqrt5 = pi * sqrt(5.0)
invSqrt6  = 2.0 * pi / sqrt(6.0)
invSqrt12 = 2.0 * pi / sqrt(12.0)

# Convert n_i to r_i using usual ehat basis
tag = [n_x - n_y,                 n_x + n_y - 2 * n_z,
       n_x + n_y + n_z - 3 * n_t, n_x + n_y + n_z + n_t]
print("tag: %d, %d, %d, %d" % (tag[0], tag[1], tag[2], tag[3]))
r = [tag[0] / sqrt(2.0),  tag[1] / sqrt(6.0),
     tag[2] / sqrt(12.0), tag[3] / sqrt(20.0)]
mag = sqrt(dot(r, r))
print("|r| = %.4g" % mag)

# Function to integrate
# Convert p_i to k_i using ghat basis
# Pattern resembles tag above, up to 2pi factors
# Integrating over dp = dk / (2pi) removes 2pi factors from measure
# Determinants now cancel, so include only factor of 2 for trace normalization
norm = 2.0 * pi * pi
def f(p):
  k = [0.0, 0.0, 0.0, 0.0]
  k[0] = (p[0] - p[1]) * invSqrt2
  k[1] = (p[0] + p[1] - 2.0 * p[2]) * invSqrt6
  k[2] = (p[0] + p[1] + p[2] - 3.0 * p[3]) * invSqrt12
  k[3] = (p[0] + p[1] + p[2] + p[3]) * halfSqrt5

  # Can only handle real part of exp
  num = exp(1.0j * (r[0] * k[0] + r[1] * k[1] + r[2] * k[2] + r[3] * k[3]))
  denom = (sin(0.5 * k[0]))**2 + (sin(0.5 * k[1]))**2 \
        + (sin(0.5 * k[2]))**2 + (sin(0.5 * k[3]))**2
  return norm * num.real / denom

integ = vegas.Integrator([[-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]])
integ(f, nitn=7, neval=Nwarm)             # Initial adaptation
result = integ(f, nitn=20, neval=Neval)   # Actual estimation

#print(result.summary())
print("result = %s    Q = %.2f" % (result, result.Q))

# Print r_I itself with propagated uncertainty
#   delta(1 / sqrt(r)) = delta(r) / 2r^(3 / 2)
# Some shenanigans necessary to get numerical results out of vegas class
temp = str(result).split('(')
val = 1.0 / sqrt(float(temp[0]))
err = float(((temp[1]).split(')'))[0])
for i in range(len(temp[0]) - 2):   # Restore leading zeroes to error
  err *= 0.1
print("--> r_I = %.8g %.4g" % (val, 0.5 * err * val**3))

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------

