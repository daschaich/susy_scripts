#!/usr/bin/python
import sys
# ------------------------------------------------------------------
# Count how many columns we can hope to get through
# from the given starting point

# Parse argument: the column at which we're starting
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<volume> <MV/sec> <start>"
  sys.exit(1)
volume = int(sys.argv[1])
matvec_per_sec = int(sys.argv[2])
iter = int(sys.argv[3])

# Total number of columns is volume * Ndat, Ndat=16DIMF
size = volume * 16 * 4    # U(2)
#size = volume * 16 * 9    # U(3)

N = 0
while N < matvec_per_sec * 24 * 3600 and iter < size:   # Total for 24 hours
 iter += 2
 toAdd = size - iter + 1
 N += toAdd
 print iter, N
if iter >= size:
  time = N / (float)(matvec_per_sec * 3600.0)
  print "Estimated time to completion:",
  print "%.1f hours" % time
# ------------------------------------------------------------------
