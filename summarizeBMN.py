#!/usr/bin/python
import os
import sys
import glob
import math
import subprocess
# ------------------------------------------------------------------
# Print out summary of all ensembles in current directory

# Make (sorted) ensemble list
ensembles = sorted(glob.glob('g*'))
if len(ensembles) < 1:
  print("No ensembles found")
  sys.exit(0)

# Print header with the following hard-coded column widths
# (not including surrounding '| ' + ' |'
#   Ensemble tag: 20 char
#   MDTU:          4 char
#   Errors:        6 char
#   Missing:       7 char
#   Therm cut:     9 char
#   Blocks:        6 char
#   poly_mod:     21 char
#   suscept:      23 char
# Total including 2x8 spaces and 7 intermediate '|': 119 char
print '|' + '{:->119}'.format('') + '|'
print '| ' + '{:20}'.format('Ensemble') + ' | MDTU',
print '| Errors | Missing | Therm cut | Blocks |',
print '{:^21}'.format('|PL|') + ' | {:^23}'.format('suscept') + ' |'
print '|' + '{:->119}'.format('') + '|'

# Cycle through and print each ensemble
for i in ensembles:
  # First check that length file exists, skip if not
  TUfile = i + '/data/TU.csv'
  if not os.path.isfile(TUfile):   # Check presence of error file
    continue

  # Now (not so efficiently) get length from last line of data/TU.csv
  for line in open(TUfile):
    temp = line
  if temp.startswith('t'):    # Zero-length ensemble
    continue
  else:                       # This ensemble is good to print
    MDTU = int(temp.split(',')[1])
    print '| {:20}'.format(i) + ' | {:>4}'.format(MDTU) + ' |',

  # Also inefficiently figure out numbers of errors and missing meas
  err = 0
  errFile = i + '/ERRORS'
  for line in open(errFile):
    err += 1
  miss = 0
  missFile = i + '/MISSING'
  for line in open(missFile):
    miss += 1
  print '{:>6}'.format(err) + ' | {:>7}'.format(miss) + ' |',

  # Also inefficiently scan therm.sh to extract therm cut
  toFind = ' ' + i + ' '        # Force exact match
  for line in open('therm.sh'):
    if i in line:
      cut = int(line.split()[2])
      break
  if cut > 0:
    print '{:>9}'.format(cut) + ' |',
  else:
    print '{:>9}'.format('unset') + ' |',

  # Get blocks and |PL| from results/poly_mod.dat if it exists
  blocks = -1
  PLfile = i + '/results/poly_mod.dat'
  if not os.path.isfile(PLfile):   # We're done
    print '{:^6}'.format('---') + ' | {:^21}'.format('---') + ' |',
    print '{:^23}'.format('---') + ' |'
  else:
    for line in open(PLfile):             # Format: dat err # blocks
      temp = line.split()
      blocks = int(temp[-1])
      dat = temp[0]
      err = temp[1]
      print '{:>6}'.format(blocks) + ' |',
      print dat + '(' + err + ') |',          # TODO: Fancify this...

  # Finally get susceptibility from results/poly_mod.suscept if it exists
  susFile = i + '/results/poly_mod.suscept'
  if os.path.isfile(susFile):   # We're done
    for line in open(susFile):            # Format: dat err # blocks
      temp = line.split()
      if not blocks == int(temp[-1]):     # Check for consistent analyses
        print "\nError: %d blocks for |PL| but %d for susceptibility" \
              % (blocks, int(temp[-1]))

      dat = temp[0]
      err = temp[1]
      print dat + '(' + err + ') |'           # TODO: Fancify this...

# Done
print '|' + '{:->119}'.format('') + '|'
# ------------------------------------------------------------------
