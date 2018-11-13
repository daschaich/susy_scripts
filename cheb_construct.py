#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Construct Chebyshev approximation to the mode number
# from the coefficients in the given output file

# Parse argument: the file to analyze
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file>")
  sys.exit(1)
toOpen = str(sys.argv[1])

# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print("ERROR:", toOpen, "does not exist")
  sys.exit(1)

# Read in coefficients
for line in open(toOpen):
  # Read in how many coefficients there are and allocate arrays
  if line.startswith('cheb_order '):
    Nterms = int((line.split())[1])
    coeffs = np.empty(Nterms, dtype = np.float)
    err = np.empty_like(coeffs)

  # Read in range of lambda, needed for variable conversion
  elif line.startswith('lambda_min '):
    la_min = float((line.split())[1])
  elif line.startswith('lambda_max '):
    la_max = float((line.split())[1])

  # Format: CHEBYSHEV c[k] coeff err
  elif line.startswith('CHEBYSHEV'):
    # First deconstruct line to find out which coeff this is
    temp = (line.split('['))[-1]        # Everything after '['
    index = int((temp.split(']'))[0])   # Number before ']'
    coeffs[index] = float((line.split())[2])
    err[index] = float((line.split())[3])

# TODO
nu = 0
numPts = 10000
la_sum = la_max + la_min
la_dif = la_max - la_min
delta = la_max / float(numPts)
for la in np.linspace(delta, la_max - delta, numPts):
  # Rescale D.Ddag eigenvalue to 'x' in range [-1, 1]
  # Recall la_max and la_min are also for D.Ddag eigenvalues |la|^2
  # Eq. 16 of arXiv:1610.01576
  x = (2.0 * la - la_sum) / la_dif
  norm = 1.0 / (np.pi * np.sqrt(1.0 - x * x))

  # rho(x) = 1 / (pi * sqrt{1 - x^2}) sum_k (2 - delta_{k0}) c_k T_k(x)
  #        = 1 / (pi * sqrt{1 - x^2}) * (2sum_k c_k T_k(x) - c_0 T_0(x))
  # Eq. 19 of arXiv:1610.01576
  rho = 2.0 * np.polynomial.chebyshev.chebval(x, coeffs)
  rho -= np.polynomial.chebyshev.chebval(x, coeffs[0])

  # TODO: Compute directly from Chebyshev polynomial rather than summing rho
  # With numPts<infty, probably missing many poles
  # Should expect severe underestimate of total 16N^2 * vol...
  rho *= norm
  nu += rho * delta
  print("%.5g %.6g %.6g" % (la, rho, nu))
# ------------------------------------------------------------------
