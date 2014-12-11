#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Compare pfaffian results from single-shot and checkpointed measurements

# Parse argument: the measurements to check
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<cfg>"
  sys.exit(1)
toCheck = int(sys.argv[1])

# First make sure we're calling this from the right place
singleShot = 'Out/check_phase.' + str(toCheck)
if not os.path.isfile(singleShot):
  print "ERROR: %s does not exist" % singleShot
  sys.exit(1)

temp = 'Out/phase_*.' + str(toCheck)
checkpointed = glob.glob(temp)
if len(checkpointed) == 0:
  print "ERROR: no files %s found" % temp
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First get and check pfaffian from single-shot measurement
for line in open(singleShot):
  # Format: Q has # columns --> # matvecs and # MBytes per core...
  if line.startswith('Q has '):
    temp = line.split()
    size = int(temp[2])   # Total number of columns
    break                 # Let's not go through the whole file just yet

# The total size of our diagonal vectors is half the number of columns
total = int(size / 2)
single_re = np.empty(total, dtype = np.float)
single_im = np.empty(total, dtype = np.float)
ckpted_re = np.empty(total, dtype = np.float)
ckpted_im = np.empty(total, dtype = np.float)

# Fill in single-shot vector
log_mag = 0.0
phase = 0.0
check_log_mag = -1.0
check_phase = -1.0
check = -1
for line in open(singleShot):
  # Format: Columns #--# of #: # matvecs in # seconds re im
  if line.startswith('Columns '):
    temp = line.split()
    i = int(int((temp[1].split('--'))[1]) / 2)  # Half of second number
    re = float(temp[9])
    im = float(temp[10])

    single_re[i - 1] = re
    single_im[i - 1] = im
    log_mag += np.log(re**2 + im**2) / 2.0
    phase += np.angle(re + 1j * im)

  # Save printed result to check
  # Format: PFAFF log_mag phase |cos| |sin|
  elif line.startswith('PFAFF '):
    temp = line.split()
    check_log_mag = float(temp[1])
    check_phase = float(temp[2])
  elif line.startswith('RUNNING COMPLETED'):
    check = 1
if check == -1:
  print singleshot, "did not complete"
  sys.exit(1)

if check_log_mag > 0:
  log_mag *= -1.0   # pf(M) = (det Q)^{-1}
  tr = np.fmod(-1.0 * phase, 2.0 * np.pi)
  if tr < 0:
    phase = tr + 2.0 * np.pi
  else:
    phase = tr
  print "Single shot roundoff check: ",
  diff = np.fabs(check_log_mag - log_mag)
  print "|%.4g - %.4g| = %.4g and" % (check_log_mag, log_mag, diff),
  diff = np.fabs(check_phase - phase)
  print "|%.4g - %.4g| = %.4g" % (check_phase, phase, diff)
else:
  print "No PFAFF line in single-shot output"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now get and check pfaffian from checkpointed measurement
log_mag = 0.0
phase = 0.0
check_log_mag = -1.0
check_phase = -1.0
for filename in checkpointed:
  for line in open(filename):
    # Format: Columns #--# of #: # matvecs in # seconds re im
    if line.startswith('Columns '):
      temp = line.split()
      i = int(int((temp[1].split('--'))[1]) / 2)  # Half of second number
      re = float(temp[9])
      im = float(temp[10])

      ckpted_re[i - 1] = re
      ckpted_im[i - 1] = im
      log_mag += np.log(re**2 + im**2) / 2.0
      phase += np.angle(re + 1j * im)

    # Save printed result to check
    # Format: PFAFF log_mag phase |cos| |sin|
    elif line.startswith('PFAFF '):
      temp = line.split()
      check_log_mag = float(temp[1])
      check_phase = float(temp[2])

if check_log_mag > 0:
  log_mag *= -1.0   # pf(M) = (det Q)^{-1}
  tr = np.fmod(-1.0 * phase, 2.0 * np.pi)
  if tr < 0:
    phase = tr + 2.0 * np.pi
  else:
    phase = tr

  print "Checkpointed roundoff check:",
  diff = np.fabs(check_log_mag - log_mag)
  print "|%.4g - %.4g| = %.4g and" % (check_log_mag, log_mag, diff),
  diff = np.fabs(check_phase - phase)
  print "|%.4g - %.4g| = %.4g" % (check_phase, phase, diff)
else:
  print "No PFAFF line in checkpointed output"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Finally compare the two element-by-element
TOL = 1e-8
count = 0
max_diff = 0.0
for i in range(total):
  re_diff = np.abs(single_re[i] - ckpted_re[i])
  im_diff = np.abs(single_im[i] - ckpted_im[i])
  if re_diff > TOL or im_diff > TOL:
    count += 1
    print "Columns %i--%i: real diff=%.4g, imag diff=%.4g" \
          % (2 * i + 1, 2 * i + 2, re_diff, im_diff)
  if re_diff > max_diff:
    max_diff = re_diff
  if im_diff > max_diff:
    max_diff = im_diff

print "\n%d of %d columns differed by more than %.1g" % (count, size, TOL)
print "The largest difference was %.4g" % max_diff
# ------------------------------------------------------------------
