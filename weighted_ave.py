#!/usr/bin/python
import os
import sys
import numpy as np
# ------------------------------------------------------------------
# Compute weighted average of data in file
# File format: dat err weight
# Also print systematic uncertainty from corresponding weighted histogram
# Parse argument: the file to analyze
if len(sys.argv) < 1:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

# Read, parse and average data
started = -1
estimates = [[], []]
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  N = len(temp) - 1
  if started < 0:
    dat = np.zeros(N, dtype = np.float)
    started = 1

  weight = float(temp[-1])
  estimates[0].append(float(temp[0]))
  estimates[1].append(weight)
  for i in range(N):
    dat[i] += float(temp[i]) * weight

if started < 0:
  print "ERROR: Didn't find any data to average"

tot_weight = sum(estimates[1])
ave = dat[0] / tot_weight
stat = dat[1] / tot_weight
print "%.6g %.4g" % (ave, stat),

# Estimate systematic uncertainty
# Step by smaller of statistical uncertainty or 0.1% of average
if stat > (ave / 1000.0):
  step = ave / 1000.0
else:
  step = stat

syst = 0.0
ratio = 0.0
while ratio < 0.68:
  syst += step
  lower = ave - syst
  upper = ave + syst
  ratio = 0.0
  for i in range(len(estimates[0])):
    if estimates[0][i] >= lower and estimates[0][i] <= upper:
      ratio += estimates[1][i]
  ratio /= tot_weight

print "%.4g" % syst
# ------------------------------------------------------------------
