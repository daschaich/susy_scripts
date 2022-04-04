#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse bosonic BMN output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting
# Normalize the Polyakov loop data by Nc and the bosonic action by Nc^2

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

ERRFILE = open('ERRORS', 'w')
MISSINGFILES = open('MISSING', 'w')

# Physical observables
SB = open('data/SB.csv', 'w')
print >> SB, "MDTU,S_B"
MYERS = open('data/Myers.csv', 'w')
print >> MYERS, "MDTU,Myers/Nt"
RATIO = open('data/ratio.csv', 'w')
print >> RATIO, "MDTU,Ratio(SO6/SO3)"
POLY = open('data/poly.csv', 'w')
print >> POLY, "ReTr(L),ImTr(L)"
POLY_MOD = open('data/poly_mod.csv', 'w')
print >> POLY_MOD, "MDTU,|Tr(L)|,ReTr(L),ImTr(L)"
SCALAR_SQUARES = open('data/scalarsquares.csv', 'w')
print >> SCALAR_SQUARES, "MDTU,Tr(X1^2),Tr(X2^2),Tr(X3^2),Tr(X4^2),Tr(X5^2),Tr(X6^2),Tr(X7^2),Tr(X8^2),Tr(X9^2)"
SCALAR_EIG_AVE = open('data/scalar_eig_ave.csv', 'w')
print >> SCALAR_EIG_AVE, "MDTU,min_ave,max_ave"
SCALAR_EIG = open('data/scalar_eig.csv', 'w')
print >> SCALAR_EIG, "MDTU,min_min,min_max,max_min,max_max"
SCALAR_EIG_WIDTHS = open('data/scalar_eig_widths.csv', 'w')
print >> SCALAR_EIG_WIDTHS, "MDTU,min_width,max_width"

# Evolution observables
ACCP = open('data/accP.csv', 'w')
print >> ACCP, "t,accP"
EXP_DS = open('data/exp_dS.csv', 'w')
print >> EXP_DS, "t,e^(-dS)"
DELTAS = open('data/deltaS.csv', 'w')
print >> DELTAS, "t,deltaS"
ABS_DS = open('data/abs_dS.csv', 'w')
print >> ABS_DS, "t,|deltaS|"
FORCE = open('data/force.csv', 'w')
print >> FORCE, "t,B"
WALLTIME = open('data/walltime.csv', 'w')
print >> WALLTIME, "t,walltime"
WALLTU = open('data/wallTU.csv', 'w')
print >> WALLTU, "t,cost"

# Run parameters
NSTEP = open('data/Nstep.csv', 'w')
print >> NSTEP, "t,N_f,N_g"
STEPSIZE = open('data/stepsize.csv', 'w')
print >> STEPSIZE, "t,eps_f,eps_g"
TLENGTH = open('data/tlength.csv', 'w')
print >> TLENGTH, "t,L"
KEY = open('data/key.csv', 'w')
print >> KEY, "t,file"
TU = open('data/TU.csv', 'w')
print >> TU, "t,MDTU"

# Status checks and running sums for the ensemble as a whole
bAct = [-1.0, -1.0]
Myers = [-1.0, -1.0]
ratio = [-1.0, -1.0]
oldcfg = 0
oldstamp = "start"
traj = 0
MDTU = 0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through files out.$load-$save from list.txt
for temp_tag in open('list.txt'):
  tag = temp_tag.rstrip()
  # Check to make sure that actual files are present
  if "Configs" in tag:
    print "No Out/out.* files found"
    print >> ERRFILE, "No Out/out.* files found"
    sys.exit(1)
  load, cfg = tag.split('-')

  # Initialize running sums and set dummy walltime
  # If the walltime isn't overwritten, then the run died
  # or its output file is corrupted
  walltime = -1
  stamp = "start"

  # Open file
  # If not found, move on to next file instead of killing whole program,
  # but print error message so I know there is a problem
  infile = 'Out/out.' + tag
  if not os.path.isfile(infile):
    print "Problem opening", infile
    print >> ERRFILE, "Problem opening", infile
    continue    # Skip to next file
#  print infile  # Monitor running status

  # If not starting from first file in this ensemble,
  # or if we seem to have skipped a file,
  # guess approximate starting trajectory
  traj_per_file = -1
  Nt = -1
  for line in open(infile):
    if line.startswith('PLACEHOLDER'):
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      walltime = -2

    # Extract Nc for bosonic action and Polyakov loop normalizations
    # Define DIMF = Nc**2 - 1
    elif line.startswith('BMN, '):
      temp = line.split(',')
      Nc = float(((temp[1]).split())[2])
      DIMF = Nc**2 - 1.0

    # Extract Nt
    elif line.startswith('nt '):
      Nt = float((line.split())[1])

    elif line.startswith('trajecs '):
      traj_per_file = int((line.split())[1])
      endtraj = traj + traj_per_file
      break       # Don't go through whole file yet

  if traj_per_file < 0:
    print infile, "never defines number of trajectories"
    print >> ERRFILE, infile, "never defines number of trajectories"
    continue    # Skip to next file
  elif Nt < 0:
    print infile, "never defines Nt"
    print >> ERRFILE, infile, "never defines Nt"
    continue    # Skip to next file

  if (traj == 0 and int(load) > 0) or (int(load) != oldcfg):
    print infile, "misses MDTU before", load
    traj = int(load)
    endtraj = traj + traj_per_file
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Cycle through lines in the file
  goodtogo = True
  for line in open(infile):
    # Check for premature termination (e.g., layout problem)
    if line.startswith('termination'):
      print infile, "reports premature termination"
      print >> ERRFILE, infile, "reports premature termination"
      goodtogo = False
      break

  if not goodtogo:
    continue        # Skip this file, which appears malformed

  # At this point we should be able to begin
  oldcfg = int(cfg)
  SCALAR_SQUARES_READY = -1     # To skip duplicate
  scalar_eig_ave = ''
  scalar_eig_ext = ''
  scalar_eig_width = ''
  for line in open(infile):
    # Extract constant run parameters
    if line.startswith('traj_length '):
      tlength = float((line.split())[1])
      print >> TLENGTH, "%d,%g" % (endtraj, tlength)
    elif line.startswith('nstep '):
      Nstep = int((line.split())[1])
      stepsize = tlength / float(Nstep)
    elif line.startswith('nstep_gauge '):
      Nstep_gauge = int((line.split())[1])
      stepsize_gauge = stepsize / float(2.0 * Nstep_gauge)
      Nstep_gauge *= 2 * Nstep
      print >> NSTEP, "%d,%d,%d" % (endtraj, Nstep, Nstep_gauge)
      print >> STEPSIZE, "%d,%g,%g" % (endtraj, stepsize, stepsize_gauge)

    elif line.startswith('Machine '):
      cpus = int((line.split())[5])

    # Check that the file loaded the appropriate configuration
    elif line.startswith('Time stamp '):
      if stamp == "start":    # Loading configuration
        stamp = line.rstrip()
        if stamp != oldstamp and oldstamp != "start":
          print infile, "time stamp doesn't match final", oldstamp
          print >> ERRFILE, infile, "time stamp doesn't match final", oldstamp
      else:                   # Saving configuration
        oldstamp = line.rstrip()    # Prepare for next file
        stamp = "start"
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Now extract evolution observables and physical observables
    # Acceptance comes before (most) measurements
    # The exception is the action
    # This is always measured twice (before and after the trajectory)
    # and which value we want depends on the accept/reject step...
    # Normalize using Nt and Nc extracted above
    # Recall DIMF = Nc**2 - 1 for SU(N)
    elif line.startswith('action: so3 '):
      if bAct[0] < 0:
        split = line.split()
        bAct[0] = float(split[10]) / (DIMF * Nt)
        Myers[0] = float(split[8]) / Nt
        ratio[0] = float(split[4]) / float(split[2])
      elif bAct[1] < 0:
        split = line.split()
        bAct[1] = float(split[10]) / (DIMF * Nt)
        Myers[1] = float(split[8]) / Nt
        ratio[1] = float(split[4]) / float(split[2])
      else:
        print infile, "lists too many action computations"
        print >> ERRFILE, infile, "lists too many action computations"

    elif ' delta S = ' in line:
      SCALAR_SQUARES_READY = 1
      traj += 1
      temp = MDTU + tlength
      MDTU = round(temp, 3)   # Round off digits in MDTU
      print >> TU, "%d,%g" % (traj, MDTU)
      print >> KEY, "%g,%s" % (traj, tag)

      dS = (line.split())[4]
      # Strip nasty comma from dS, if necessary
      if ',' in dS:
        dS = (dS.split(','))[0]
      print >> DELTAS, "%d,%g" % (traj, float(dS))
      print >> EXP_DS, "%d,%g" % (traj, math.exp(-1.0 * float(dS)))

      # For RMS, don't have an easy way to average over many measurements
      # Instead just print out absolute value and consider its running average
      print >> ABS_DS, "%d,%g" % (traj, abs(float(dS)))

      # Acceptance is smeared out by running averages
      if line.startswith('ACCEPT'):
        print >> ACCP, "%d,1" % traj
        print >> SB, "%g,%g" % (MDTU, bAct[1])
        print >> MYERS, "%g,%g" % (MDTU, Myers[1])
        print >> RATIO, "%g,%g" % (MDTU, ratio[1])
      else:
        print >> ACCP, "%d,0" % traj
        print >> SB, "%g,%g" % (MDTU, bAct[0])
        print >> MYERS, "%g,%g" % (MDTU, Myers[0])
        print >> RATIO, "%g,%g" % (MDTU, ratio[0])
      bAct = [-1.0, -1.0]                          # Reset
      Myers = [-1.0, -1.0]
      ratio = [-1.0, -1.0]

    # Forces -- take maxima rather than average if possible
    elif line.startswith('MONITOR_FORCE_GAUGE '):
      force_g = float((line.split())[-1])
      print >> FORCE, "%d,%g" % (traj, force_g)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Polyakov loop
    # First normalized using Nc extracted above
    elif line.startswith('GMES '):
      temp = line.split()
      poly_r = float(temp[1]) / Nc
      poly_i = float(temp[2]) / Nc
      print >> POLY, "%g,%g" % (poly_r, poly_i)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD, "%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Nine scalar squares
    # Skip duplicate measurement before first trajectory
    elif line.startswith('SCALAR SQUARES ') and SCALAR_SQUARES_READY > 0:
      temp = line.split()
      X1 = float(temp[2])
      X2 = float(temp[3])
      X3 = float(temp[4])
      X4 = float(temp[5])
      X5 = float(temp[6])
      X6 = float(temp[7])
      X7 = float(temp[8])
      X8 = float(temp[9])
      X9 = float(temp[10])
      print >> SCALAR_SQUARES, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g" \
                               % (MDTU, X1, X2, X3, X4, X5, X6, X7, X8, X9)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Scalar eigenvalues
    # Only look at largest and smallest (most negative) to keep this manageable
    elif line.startswith('SCALAR_EIG '):
      temp = line.split()
      index = int(temp[1])
      if index == 0 or index == Nc - 1:
        scalar_eig_ave += ',' + str(temp[2])
        scalar_eig_ext += ',' + str(temp[4]) + ',' + str(temp[5])
        scalar_eig_width += ',' + str(temp[3])
      if index == Nc - 1:
        print >> SCALAR_EIG_AVE, "%g%s" % (MDTU, scalar_eig_ave)
        print >> SCALAR_EIG, "%g%s" % (MDTU, scalar_eig_ext)
        print >> SCALAR_EIG_WIDTHS, "%g%s" % (MDTU, scalar_eig_width)
        scalar_eig_ave = ''
        scalar_eig_ext = ''
        scalar_eig_width = ''
    # ------------------------------------------------------------

    # Store total walltime to average at the end
    elif line.startswith('Time = '):
      walltime = float((line.split())[2])
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Check to see if run seems to have finished properly
  if walltime == -1:
    print infile, "didn't print final timing"
    print >> ERRFILE, infile, "didn't print final timing"
  elif walltime == -2:
    # Placeholder file -- error has been addressed as well as possible,
    # but don't print nonsense wall clock time
    bAct = [-1.0, -1.0]                          # Reset
    pass
  else:   # We are good to go
    ave_time = walltime / traj_per_file
    print >> WALLTIME, "%d,%g" % (traj, ave_time)
    TUtime = walltime * float(cpus) / (60.0 * tlength * traj_per_file)
    print >> WALLTU, "%d,%g" % (traj, TUtime)
  # ----------------------------------------------------------------
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ERRFILE.close()
MISSINGFILES.close()
SB.close()
MYERS.close()
RATIO.close()
POLY.close()
POLY_MOD.close()
SCALAR_SQUARES.close()
SCALAR_EIG_AVE.close()
SCALAR_EIG.close()
SCALAR_EIG_WIDTHS.close()
ACCP.close()
EXP_DS.close()
DELTAS.close()
ABS_DS.close()
FORCE.close()
WALLTIME.close()
WALLTU.close()
NSTEP.close()
STEPSIZE.close()
TLENGTH.close()
KEY.close()
TU.close()
# ------------------------------------------------------------------
