#!/usr/bin/python
from __future__ import print_function     # For vegas
import sys
import time
from numpy import sqrt, sin, cos, pi
import vegas
# ------------------------------------------------------------------
# Print r_lat (as defined in notes)
# corresponding to input integer displacement three-vector

# Parse arguments: Three-component integer vector
if len(sys.argv) < 5:
  print("Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <Neval>")
  sys.exit(1)
n_x = int(sys.argv[1])
n_y = int(sys.argv[2])
n_z = int(sys.argv[3])
Neval = int(sys.argv[4])
Nwarm = int(Neval / 10)
runtime = -time.time()

# Convert (n_x, n_y, n_z) to (r_1, r_2, r_3)
r = [(n_x - n_y) / sqrt(2.0),           \
     (n_x + n_y - 2 * n_z) / sqrt(6.0), \
     (n_x + n_y + n_z) / sqrt(12.0)]
mag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
print("|r| = %.4g" % mag)

# Function for n = (1, 1, 1)
def f(k):
  num = cos(r[0] * k[0]) * cos(r[1] * k[1]) * cos(r[2] * k[2])
  denom = (sin(0.5 * k[0]))**2 + (sin(0.5 * k[1]))**2 + (sin(0.5 * k[2]))**2
  return num / (8.0 * pi**2 * denom)

integ = vegas.Integrator([[-1.0 * pi, pi], [-1.0 * pi, pi], [-1.0 * pi, pi]])
integ(f, nitn=7, neval=Nwarm)             # Initial adaptation
result = integ(f, nitn=20, neval=Neval)   # Actual estimation

#print(result.summary())
print("result = %s    Q = %.2f" % (result, result.Q))

# Print r_lat itself with propagated uncertainty delta(1 / r) = delta(r) / r^2
# Some shenanigans necessary to get numerical results out of vegas class
temp = str(result).split('(')
val = 1.0 / float(temp[0])
err = float(((temp[1]).split(')'))[0])
for i in range(len(temp[0]) - 2):   # Restore leading zeroes to error
  err *= 0.1
print("--> r_lat = %.8g %.4g" % (val, err * val * val))

runtime += time.time()
print("Runtime: %0.1f seconds" % runtime)
# ------------------------------------------------------------------

