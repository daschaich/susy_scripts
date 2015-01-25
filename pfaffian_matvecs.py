#!/usr/bin/python
import sys
# ------------------------------------------------------------------
# Count how many columns we can hope to get through
# from the given starting point

# Parse argument: the column at which we're starting
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<N> <volume> <MV/sec> <start>"
  sys.exit(1)
Nc = int(sys.argv[1])
volume = int(sys.argv[2])
matvec_per_sec = int(sys.argv[3])
iter = int(sys.argv[4])

# Total number of columns is volume * Ndat, Ndat=16DIMF
size = volume * 16 * Nc * Nc

N = 0
while N < matvec_per_sec * 24 * 3600 and iter < size:   # Total for 24 hours
 iter += 2
 toAdd = size - iter + 1
 N += toAdd
if iter >= size:
  time = N / (float)(matvec_per_sec * 3600.0)
  print "Estimated time to completion:",
  print "%.1f hours" % time
else:
  print "Should be able to reach %d (%d matvecs) in 24 hours" \
        % (iter, N)
# ------------------------------------------------------------------
