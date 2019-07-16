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

# Print header with the following eight hard-coded column widths
# for the ensemble tag, MDTU, errors, missing, therm cut,
#         blocks, |PL| and susceptibility
# Total includes 2x8 spaces and 7 intermediate '|': 109 char
widths = [20, 4, 6, 7, 9, 6, 11, 12]
tot_width = sum(widths) + 2 * 8 + 7
print '|' + '{:->{w}}'.format('', w=tot_width) + '|'
print '| ' + '{:{w}}'.format('Ensemble', w=widths[0]) + ' | MDTU',
print '| Errors | Missing | Therm cut | Blocks |',
print '{:^{w}}'.format('|PL|', w=widths[6]) + ' |',
print '{:^{w}}'.format('suscept', w=widths[7]) + ' |'
print '|' + '{:->{w}}'.format('', w=tot_width) + '|'

# Convenience quantities for |PL| and susceptibility formatting
PLwidth = 5
PLscale = 1e5
susWidth = 6
susScale = 1e6

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
    print '| {:{w}}'.format(i, w=widths[0]) + ' |',
    print '{:>{w}}'.format(MDTU, w=widths[1]) + ' |',

  # Also inefficiently figure out numbers of errors and missing meas
  err = 0
  errFile = i + '/ERRORS'
  for line in open(errFile):
    err += 1
  miss = 0
  missFile = i + '/MISSING'
  for line in open(missFile):
    miss += 1
  print '{:>{w}}'.format(err, w=widths[2]) + ' |',
  print '{:>{w}}'.format(miss, w=widths[3]) + ' |',

  # Also inefficiently scan therm.sh to extract therm cut
  toFind = ' ' + i + ' '        # Force exact match
  for line in open('therm.sh'):
    if i in line:
      cut = int(line.split()[2])
      break
  if cut > 0:
    print '{:>{w}}'.format(cut, w=widths[4]) + ' |',
  else:
    print '{:>{w}}'.format('unset', w=widths[4]) + ' |',

  # Get blocks and |PL| from results/poly_mod.dat if it exists
  blocks = -1
  PLfile = i + '/results/poly_mod.dat'
  if not os.path.isfile(PLfile):   # We're done
    print '{:^{w}}'.format('---', w=widths[5]) + ' |',
    print '{:^{w}}'.format('---', w=widths[6]) + ' |',
    print '{:^{w}}'.format('---', w=widths[7]) + ' |'
  else:
    for line in open(PLfile):             # Format: dat err # blocks
      temp = line.split()
      blocks = int(temp[-1])
      dat = float(temp[0])
      err = int(PLscale * float(temp[1]))
      print '{:>{w}}'.format(blocks, w=widths[5]) + ' |',
      print '{:>.{w}f}'.format(dat, w=PLwidth) + '(' + str(err) + ') |',

  # Finally get susceptibility from results/poly_mod.suscept if it exists
  susFile = i + '/results/poly_mod.suscept'
  if os.path.isfile(susFile):   # We're done
    for line in open(susFile):            # Format: dat err # blocks
      temp = line.split()
      if not blocks == int(temp[-1]):     # Check for consistent analyses
        print "\nError: %d blocks for |PL| but %d for susceptibility" \
              % (blocks, int(temp[-1]))

      dat = float(temp[0])
      err = int(susScale * float(temp[1]))
      print '{:>.{w}f}'.format(dat, w=susWidth) + '(' + str(err) + ') |'

# Done
print '|' + '{:->{w}}'.format('', w=tot_width) + '|'
# ------------------------------------------------------------------
