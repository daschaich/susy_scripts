#!/usr/bin/python
from __future__ import print_function     # For vegas
import sys
import time
import numpy as np
import vegas
# ------------------------------------------------------------------
# Print r_I (as defined in notes)
# corresponding to input integer displacement four-vector
# Use vegas to numerically calculate the infinite-volume Fourier transform

# Try to set up parallel evaluation
# Need some global variables (initialized in main()
n = [0, 0, 0, 0]
twopi = 2.0 * np.pi
piSq = np.pi * np.pi
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Function to integrate -- can only handle real part of exp
# Also can't take the square root since negative evaluations seem possible
# Don't add any extra normalization (determinants cancel)
# Need to handle multiple sets of p's at once...
def f(p):
  num = np.zeros(p.shape[0], dtype = np.float)
  denom = np.zeros_like(num)
  temp = np.zeros_like(num)
  for i in range(p.shape[0]):
    temp[i] = n[0] * p[i][0] + n[1] * p[i][1] \
            + n[2] * p[i][2] + n[3] * p[i][3]
    num[i] = np.cos(twopi * temp[i])
    denom[i] = (np.sin(np.pi * p[i][0]))**2 + (np.sin(np.pi * p[i][1]))**2 \
             + (np.sin(np.pi * p[i][2]))**2 + (np.sin(np.pi * p[i][3]))**2
  return piSq * num / denom
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: Four-component integer vector
# and number of evaluations in vegas calculation
def main():
  if len(sys.argv) < 6:
    print("Usage:", str(sys.argv[0]), "<n_x> <n_y> <n_z> <n_t> <Neval>")
    sys.exit(1)
  n[0] = int(sys.argv[1])
  n[1] = int(sys.argv[2])
  n[2] = int(sys.argv[3])
  n[3] = int(sys.argv[4])
  Neval = int(sys.argv[5])
  Nwarm = int(Neval / 10)
  runtime = -time.time()

  # A4* basis vectors as 4x4 matrix
  invSqrt2  = 1.0 / np.sqrt(2.0)
  invSqrt6  = 1.0 / np.sqrt(6.0)
  invSqrt12 = 1.0 / np.sqrt(12.0)
  invSqrt20 = 1.0 / np.sqrt(20.0)
  ehat = [[invSqrt2,        invSqrt6,        invSqrt12,        invSqrt20], \
          [-1.0 * invSqrt2, invSqrt6,        invSqrt12,        invSqrt20], \
          [0,               -2.0 * invSqrt6, invSqrt12,        invSqrt20], \
          [0,               0,               -3.0 * invSqrt12, invSqrt20]]

  # Compute naive |r| using ehat (print after defining fparallel)
  r = [0.0, 0.0, 0.0, 0.0]
  for i in range(4):
    for j in range(4):
      r[j] += n[i] * ehat[i][j]
  mag = np.sqrt(np.dot(r, r))

  # Manage integration
  # Integrating over dp = dk / (2pi) removes 2pi factors from measure
  integ = vegas.Integrator(4 * [[0.0, 1.0]])
  # Convert to MPI integrand
  fparallel = vegas.MPIintegrand(f)
  if fparallel.rank == 0:
    print("|r| = %.4g" % mag)
  integ(fparallel, nitn=7, neval=Nwarm)             # Initial adaptation
  result = integ(fparallel, nitn=20, neval=Neval)   # Actual estimation

  if fparallel.rank == 0:
    #print(result.summary())
    print("result = %s    Q = %.2f" % (result, result.Q))

  # Print r_I itself with propagated uncertainty
  #   delta(1 / sqrt(r)) = delta(r) / 2r^(3 / 2)
  # Some shenanigans necessary to get numerical results out of vegas class
  temp = str(result).split('(')
  val = 1.0 / np.sqrt(float(temp[0]))
  err = float(((temp[1]).split(')'))[0])
  for i in range(len(temp[0]) - 2):   # Restore leading zeroes to error
    err *= 0.1
  if fparallel.rank == 0:
    print("--> r_I = %.8g %.4g" % (val, 0.5 * err * val**3))

  runtime += time.time()
  if fparallel.rank == 0:
    print("Runtime: %0.1f seconds" % runtime)

if __name__ == '__main__':
  main()
# ------------------------------------------------------------------
