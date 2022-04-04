#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse BMN output files for a single ensemble,
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
SF = open('data/SF.csv', 'w')
print >> SF, "MDTU,S_F"
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
EIG = open('data/eig.csv', 'w')
print >> EIG, "MDTU,min_eig"

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
print >> FORCE, "t,G,F1,F2,F3"
CG_ITERS = open('data/cg_iters.csv', 'w')
print >> CG_ITERS, "t,cg_iters"
WALLTIME = open('data/walltime.csv', 'w')
print >> WALLTIME, "t,walltime"
WALLTU = open('data/wallTU.csv', 'w')
print >> WALLTU, "t,cost"
COND_NUM = open('data/cond_num.csv', 'w')
print >> COND_NUM, "t,cond_num"

# Run parameters
NSTEP = open('data/Nstep.csv', 'w')
print >> NSTEP, "t,N_f,N_g"
STEPSIZE = open('data/stepsize.csv', 'w')
print >> STEPSIZE, "t,eps_f,eps_g"
NORDER = open('data/Norder.csv', 'w')
print >> NORDER, "t,Norder"
TLENGTH = open('data/tlength.csv', 'w')
print >> TLENGTH, "t,L"
KEY = open('data/key.csv', 'w')
print >> KEY, "t,file"
TU = open('data/TU.csv', 'w')
print >> TU, "t,MDTU"

# Status checks and running sums for the ensemble as a whole
fermAct = [-1.0, -1.0]
bAct = [-1.0, -1.0]
Myers = [-1.0, -1.0]
ratio = [-1.0, -1.0]
oldcfg = 0
oldstamp = "start"
CG = 1
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
  Nroot = 1   # Default
  min_eig = 1
  max_eig = -1
  SCALAR_SQUARES_READY = -1     # To skip duplicate
  scalar_eig_ave = ''
  scalar_eig_ext = ''
  scalar_eig_width = ''
  for line in open(infile):
    # See how many fermion forces we will have below
    # Retain case insensitivity for now
    if line.lower().startswith('using nroot '):
      Nroot = int((line.split())[3])

    # Extract spectral range for eigenvalues
    # Format: RHMC Norder # for spectral range [min, max]
    elif line.startswith('RHMC Norder '):
      if 'spectral' in line:
        temp = line.rstrip()       # Kill newline
        Norder = int(temp.split()[2])
        temp2 = temp.rstrip(']')   # Kill ]
        temp = (temp2.split('['))[-1]
        temp2 = temp.split(',')
        min_eig = float(temp2[0])
        max_eig = float(temp2[1])
        print >> NORDER, "%d,%d" % (endtraj, Norder)
      else:         # Original 15-pole format didn't state spectral range
        min_eig = 1.0e-7
        max_eig = 1000.0
        print >> NORDER, "%d,%d" % (endtraj, 15)

    # Extract constant run parameters
    elif line.startswith('traj_length '):
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
    # The exception is the fermion action
    # This is always measured twice (before and after the trajectory)
    # and which value we want depends on the accept/reject step...
    # Normalize using Nt and Nc extracted above
    # Recall DIMF = Nc**2 - 1 for SU(N)
    elif line.startswith('action: so3 '):
      if fermAct[0] < 0:
        split = line.split()
        fermAct[0] = float(split[12]) / (16.0 * Nt * DIMF)
        bAct[0] = float(split[10]) / (DIMF * Nt)
        Myers[0] = float(split[8]) / Nt
        ratio[0] = float(split[4]) / float(split[2])
      elif fermAct[1] < 0:
        split = line.split()
        fermAct[1] = float(split[12]) / (16.0 * Nt * DIMF)
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
        print >> SF, "%d,%g" % (MDTU, fermAct[1])   # New action
        print >> SB, "%g,%g" % (MDTU, bAct[1])
        print >> MYERS, "%g,%g" % (MDTU, Myers[1])
        print >> RATIO, "%g,%g" % (MDTU, ratio[1])
      else:
        print >> ACCP, "%d,0" % traj
        print >> SF, "%d,%g" % (MDTU, fermAct[0])   # Original action
        print >> SB, "%g,%g" % (MDTU, bAct[0])
        print >> MYERS, "%g,%g" % (MDTU, Myers[0])
        print >> RATIO, "%g,%g" % (MDTU, ratio[0])
      fermAct = [-1.0, -1.0]                        # Reset
      bAct = [-1.0, -1.0]
      Myers = [-1.0, -1.0]
      ratio = [-1.0, -1.0]

    # Forces -- take maxima rather than average if possible
    elif line.startswith('MONITOR_FORCE_GAUGE '):
      force_g = float((line.split())[-1])
    elif line.startswith('MONITOR_FORCE_FERMION0 '):
      force_f = float((line.split())[-1])
      if Nroot == 1:
        print >> FORCE, "%d,%g,%g,null,null" % (traj, force_g, force_f)
    elif line.startswith('MONITOR_FORCE_FERMION1 '):
      force_f2 = float((line.split())[-1])
      if Nroot == 2:
        print >> FORCE, "%d,%g,%g,%g,null" % (traj, force_g, force_f, force_f2)
    elif line.startswith('MONITOR_FORCE_FERMION2 '):
      force_f3 = float((line.split())[-1])
      if Nroot == 3:
        print >> FORCE, "%d,%g,%g,%g,%g" \
                        % (traj, force_g, force_f, force_f2, force_f3)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Polyakov loop and CG iterations
    # First normalized using Nc extracted above
    elif line.startswith('GMES '):
      temp = line.split()
      print >> CG_ITERS, "%g,%g" % (traj, float(temp[3]))

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

    # Check to make sure CG always converged
    elif 'CONGRAD' in line:
      CG = -1

    # Store total walltime to average at the end
    elif line.startswith('Time = '):
      walltime = float((line.split())[2])
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Check to see if run seems to have finished properly
  if CG == -1:
    print infile, "encountered CG non-convergence"
    print >> ERRFILE, infile, "encountered CG non-convergence"
    CG = 1
  if walltime == -1:
    print infile, "didn't print final timing"
    print >> ERRFILE, infile, "didn't print final timing"
  elif walltime == -2:
    # Placeholder file -- error has been addressed as well as possible,
    # but don't print nonsense wall clock time
    fermAct = [-1.0, -1.0]                        # Reset
    pass
  else:   # We are good to go
    ave_time = walltime / traj_per_file
    print >> WALLTIME, "%d,%g" % (traj, ave_time)
    TUtime = walltime * float(cpus) / (60.0 * tlength * traj_per_file)
    print >> WALLTU, "%d,%g" % (traj, TUtime)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "eig" file, if it is present
  # These seem to come in sets of 16 (which we don't check too carefully)
  # For now only consider the smallest
  infile = 'Out/eig.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
  else:
    check = -1    # Check whether file completed successfully

    # Go
    for line in open(infile):
      if line.startswith('Time stamp '):
        stamp = line.rstrip()
        if stamp != oldstamp:
          print infile, "time stamp doesn't match final", oldstamp
          print >> ERRFILE, infile, "time stamp doesn't match final", oldstamp

      elif line.startswith('EIGENVALUE 0 '):
        eig = float((line.split())[2])
        if eig < min_eig:        # Check spectral range
          print infile, "exceeds RHMC spectral range:",
          print "%.4g not in [%.4g, %.4g]" % (eig, min_eig, max_eig)
          print >> ERRFILE, infile, "exceeds RHMC spectral range:",
          print >> ERRFILE, "%.4g not in [%.4g, %.4g]" % (eig, min_eig, max_eig)

      elif line.startswith('BIGEIGVAL  0 '):    # Check spectral range
        dat = float((line.split())[2])
        if dat > max_eig:
          print infile, "exceeds RHMC spectral range:",
          print "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)
          print >> ERRFILE, infile, "exceeds RHMC spectral range:",
          print >> ERRFILE, "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)

        # Monitor (log of) condition number
        cond_num = math.log(dat / eig)
        print >> COND_NUM, "%d,%g" % (traj, cond_num)

      elif 'WARNING' in line:
        print infile, "saturated eigenvalue iterations"
        print >> ERRFILE, infile, "saturated eigenvalue iterations"
      elif line.startswith('RUNNING COMPLETED'):
        if check == 1:    # Check that we have one measurement per file
          print infile, "reports two measurements"
          print >> ERRFILE, infile, "reports two measurements"
        check = 1
    if check == -1:
      print infile, "did not complete"
      print >> ERRFILE, infile, "did not complete"

    print >> EIG, "%g,%g" % (MDTU, eig)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ERRFILE.close()
MISSINGFILES.close()
SB.close()
MYERS.close()
RATIO.close()
SF.close()
POLY.close()
POLY_MOD.close()
SCALAR_SQUARES.close()
SCALAR_EIG_AVE.close()
SCALAR_EIG.close()
SCALAR_EIG_WIDTHS.close()
EIG.close()
ACCP.close()
EXP_DS.close()
DELTAS.close()
ABS_DS.close()
FORCE.close()
CG_ITERS.close()
WALLTIME.close()
WALLTU.close()
COND_NUM.close()
NSTEP.close()
STEPSIZE.close()
NORDER.close()
TLENGTH.close()
KEY.close()
TU.close()
# ------------------------------------------------------------------
