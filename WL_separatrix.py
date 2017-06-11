#!/usr/bin/env python
import os
import sys
import glob
# ------------------------------------------------------------------
# Compute interpolator between confined and deconfined behavior
# For now only consider unitarized (and not full) Wilson lines
# TODO: Need an estimate of statistical uncertainties...

# Parse arguments: First is thermalization cut,
# second is radius of separatrix
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<cut> <r>"
  sys.exit(1)
cut = int(sys.argv[1])
radius = float(sys.argv[2])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Extract Nc from path -- it's after 'Nc' then before '_'
cwd = os.getcwd()
temp = (cwd.split('Nc'))[1]
Nc = int((temp.split('_'))[0])

# Extract tag from path as everything after the last '/'
# Then extract rt
tag = (cwd.split('/'))[-1]
temp = (tag.split('_'))[0]
rt = (temp.split('rt'))[1]

# Extract Nc-normalized modulus as third datum on each line of data file
# (Format: MDTU,|Lx|,|Ly|,|Lz|,|L5|)
low = 0
high = 0
for line in open('data/lines_mod_polar.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  dat = float(temp[3])
  if dat < radius:
    low += 1
  else:
    high += 1
tot = low + high

# Compute S = [(Nc + 1) * low - tot] / [(Nc - 1) * low + tot]
# Print number of measurements to allow offline checks
num = (Nc + 1) * low - tot
den = (Nc - 1) * low + tot
s = float(num) / float(den)
print "%s %.4g #" % (rt, s),
print "%d + %d = %d measurements for %s" % (low, high, tot, tag)
# ------------------------------------------------------------------
