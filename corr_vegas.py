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
n = [0, 0, 0, 0]
n[0] = int(sys.argv[1])
n[1] = int(sys.argv[2])
n[2] = int(sys.argv[3])
n[3] = int(sys.argv[4])
Neval = int(sys.argv[5])
Nwarm = int(Neval / 10)
runtime = -time.time()

# A4* basis vectors as 4x4 matrix
invSqrt2  = 1.0 / sqrt(2.0)
invSqrt6  = 1.0 / sqrt(6.0)
invSqrt12 = 1.0 / sqrt(12.0)
invSqrt20 = 1.0 / sqrt(20.0)
twopi     = 2.0j * pi
piSq      = pi * pi
ehat = [[invSqrt2,        invSqrt6,        invSqrt12,        invSqrt20], \
        [-1.0 * invSqrt2, invSqrt6,        invSqrt12,        invSqrt20], \
        [0,               -2.0 * invSqrt6, invSqrt12,        invSqrt20], \
        [0,               0,               -3.0 * invSqrt12, invSqrt20]]

# Print naive |r| using ehat
r = [0.0, 0.0, 0.0, 0.0]
for i in range(4):
  for j in range(4):
    r[j] += n[i] * ehat[i][j]

mag = sqrt(dot(r, r))
print("|r| = %.4g" % mag)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Function to integrate -- can only handle real part of exp
# Also can't take the square root since negative evaluations seem possible
# Integrating over dp = dk / (2pi) removes 2pi factors from measure
# Don't add any extra normalization (determinants cancel)
def f(p):
  num = exp(twopi * (n[0] * p[0] + n[1] * p[1] + n[2] * p[2] + n[3] * p[3]))
  denom = (sin(pi * p[0]))**2 + (sin(pi * p[1]))**2 \
        + (sin(pi * p[2]))**2 + (sin(pi * p[3]))**2
  return piSq * num.real / denom

integ = vegas.Integrator(4 * [[0.0, 1.0]])
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
