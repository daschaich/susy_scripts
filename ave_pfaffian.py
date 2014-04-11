#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compute <e^{i alpha)> from pfaffian phase measurements
# not doing any cutting or blocking at the moment

# No arguments
if len(sys.argv) > 1:
  print "Usage:", str(sys.argv[0])
  sys.exit(1)

# Cycle through both single-shot and checkpointed results
# but don't waste time on intermediate checkpoints
dat = []
files = glob.glob('phase.*')
files += glob.glob('phase_*-end.*')
if len(files) == 0:
  print "ERROR: no files phase.* phase_*-end.* in this directory"
  sys.exit(1)

for filename in files:
  for line in open(filename):
    # Format: PFAFF log_mag phase |cos(phase)| |sin(phase)|
    if line.startswith('PFAFF '):
      temp = line.split()
      dat.append(np.exp(1.j * float(temp[2])))

# Construct jackknife samples through single elimination
tot = np.sum(dat)
N = len(dat)
jkre = np.empty(N, dtype = np.float)
jkim = np.empty(N, dtype = np.float)
jkmag = np.empty(N, dtype = np.float)
jkphase = np.empty(N, dtype = np.float)
for i in range(N):
  temp = (tot - dat[i]) / (N - 1.0)
  jkre[i] = np.real(temp)
  jkim[i] = np.imag(temp)
  jkmag[i] = np.absolute(temp)
  jkphase[i] = np.angle(temp)

ave_re = np.average(jkre)
var_re = (N - 1.0) * np.sum((jkre - ave_re)**2) / float(N)

ave_im = np.average(jkim)
var_im = (N - 1.0) * np.sum((jkim - ave_im)**2) / float(N)

ave_mag = np.average(jkmag)
var_mag = (N - 1.0) * np.sum((jkmag - ave_mag)**2) / float(N)

ave_phase = np.average(jkphase)
var_phase = (N - 1.0) * np.sum((jkphase - ave_phase)**2) / float(N)

print "From %d measurements," % N,
print "<e^{i alpha}> = (%.4g, %.4g) = %.4g e^{%.4gi}" % \
      (ave_re, ave_im, ave_mag, ave_phase)
print "Real: %.6g +/- %.4g;" % (ave_re, np.sqrt(var_re)),
print "Imag: %.6g +/- %.4g;" % (ave_im, np.sqrt(var_im))
print "Mag: %.6g +/- %.4g;" % (ave_mag, np.sqrt(var_mag)),
print "Phase: %.6g +/- %.4g;" % (ave_phase, np.sqrt(var_phase))
# ------------------------------------------------------------------
