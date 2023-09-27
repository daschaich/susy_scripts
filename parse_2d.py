#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse N=(2, 2) SYM output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting
# Normalize the Wilson line data by Nc and the bosonic action by 3Nc^2 / 2

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

ERRFILE = open('ERRORS', 'w')
MISSINGFILES = open('MISSING', 'w')

# Physical observables
PLAQ = open('data/plaq.csv', 'w')
print >> PLAQ, "MDTU,plaq"
SB = open('data/SB.csv', 'w')
print >> SB, "MDTU,S_B"
SF = open('data/SF.csv', 'w')
print >> SF, "MDTU,S_F"
POLY = open('data/poly.csv', 'w')
print >> POLY, "ReTr(L),ImTr(L)"
POLY_POLAR = open('data/poly_polar.csv', 'w')
print >> POLY_POLAR, "ReTr(L),ImTr(L)"
POLY_MOD = open('data/poly_mod.csv', 'w')
print >> POLY_MOD, "MDTU,|Tr(L)|,ReTr(L),ImTr(L)"
POLY_MOD_POLAR = open('data/poly_mod_polar.csv', 'w')
print >> POLY_MOD_POLAR, "MDTU,|Tr(L)|,ReTr(L),ImTr(L)"
LINES = open('data/lines.csv', 'w')
print >> LINES, "ReTr(Lx),ImTr(Lx)"
LINES_POLAR = open('data/lines_polar.csv', 'w')
print >> LINES_POLAR, "ReTr(Lx),ImTr(Lx)"
LINES_MOD = open('data/lines_mod.csv', 'w')
print >> LINES_MOD, "MDTU,|Tr(Lx)|"
LINES_MOD_POLAR = open('data/lines_mod_polar.csv', 'w')
print >> LINES_MOD_POLAR, "MDTU,|Tr(Lx)|"
FLINK = open('data/Flink.csv', 'w')
print >> FLINK, "MDTU,link"
DET = open('data/det.csv', 'w')
print >> DET, "MDTU,|det-1|^2,1-Re(det),Im(det)"
WIDTHS = open('data/widths.csv', 'w')
print >> WIDTHS, "MDTU,plaq,Re(det),Im(det),link"
SCALAR_EIG_AVE = open('data/scalar_eig_ave.csv', 'w')
print >> SCALAR_EIG_AVE, "MDTU,min_ave,max_ave"
SCALAR_EIG = open('data/scalar_eig.csv', 'w')
print >> SCALAR_EIG, "MDTU,min_min,min_max,max_min,max_max5"
SCALAR_EIG_WIDTHS = open('data/scalar_eig_widths.csv', 'w')
print >> SCALAR_EIG_WIDTHS, "MDTU,min_width,max_width"
BILIN = open('data/bilin.csv', 'w')
print >> BILIN, "MDTU,susyTrans,Im(bilin)"
EIG = open('data/eig.csv', 'w')
print >> EIG, "MDTU,0,2,4,6,8,10"
BIGEIG = open('data/bigeig.csv', 'w')
print >> BIGEIG, "MDTU,0,2,4,6,8,10"
ENERGY = open('data/energy.csv', 'w')
print >> ENERGY, "MDTU,energy"
SCALARSQ = open('data/scalar_sq.csv', 'w')
print >> SCALARSQ, "MDTU,tr_xsq"

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
  vol = -1
  lam = -1
  for line in open(infile):
    if line.startswith('PLACEHOLDER'):
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      walltime = -2

    # Extract Nc for bosonic action and Polyakov loop normalizations
    elif line.startswith('N=(2,'):
      temp = line.split(',')
      Nc = float(((temp[2]).split())[2])
      DIMF = Nc**2
    # Extract volume for fermion action normalization
    elif line.startswith('nx '):
      vol = float((line.split())[1])
    elif line.startswith('nt '):
      vol *= float((line.split())[1])

    # Make sure no overflow errors parsing seed
    elif line.startswith('iseed '):
      if int((line.split())[1]) == -1:
        print infile, "has suspicious seed -1"
        print >> ERRFILE, infile, "has suspicious seed -1"

    elif line.startswith('trajecs '):
      traj_per_file = int((line.split())[1])
      endtraj = traj + traj_per_file
    elif line.startswith('lambda '):
      lam = float((line.split())[1])
      break       # Don't go through whole file yet

  if traj_per_file < 0:
    print infile, "never defines number of trajectories"
    print >> ERRFILE, infile, "never defines number of trajectories"
    continue    # Skip to next file
  elif vol < 0:
    print infile, "never defines lattice volume"
    print >> ERRFILE, infile, "never defines lattice volume"
    continue    # Skip to next file
  elif lam < 0:
    print infile, "never defines lattice lambda"
    print >> ERRFILE, infile, "never defines lattice lambda"
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

    # Check for unitarity problem
    # Report for info, but don't skip if found
    if line.startswith('Unitarity problem'):
      print infile, "reports unitarity problem"
      print >> ERRFILE, infile, "reports unitarity problem"
      break

  if not goodtogo:
    continue        # Skip this file, which appears malformed

  # At this point we should be able to begin
  oldcfg = int(cfg)
  fresh = -1
  Nroot = 1   # Default
  min_eig = 1
  max_eig = -1
  NEED_LINES = 1    # May not need to check measurement file
  NEED_DET = 1      # for these if they're in the main output file
  NEED_WIDTHS = 1
  NEED_SCALAR_EIGS = 1
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

    # Check whether this is a fresh start
    elif line.startswith('fresh'):
      fresh = 1

    elif line.startswith('Machine '):
      cpus = int((line.split())[5])

    # Check that the file loaded the appropriate configuration
    elif line.startswith('Time stamp '):
      if stamp == "start":    # Loading configuration
        stamp = line.rstrip()
        if fresh > 0:         # Ensure checks pass for fresh start
          oldstamp = stamp
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
    # Normalize using volume and Nc extracted above
    elif line.startswith('action: gauge '):
      if fermAct[0] < 0:
        fermAct[0] = float((line.split())[8]) / (4.0 * vol * DIMF)
      elif fermAct[1] < 0:
        fermAct[1] = float((line.split())[8]) / (4.0 * vol * DIMF)
      else:
        print infile, "lists too many action computations"
        print >> ERRFILE, infile, "lists too many action computations"

    elif ' delta S = ' in line:
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
      else:
        print >> ACCP, "%d,0" % traj
        print >> SF, "%d,%g" % (MDTU, fermAct[0])   # Original action
      fermAct = [-1.0, -1.0]                        # Reset

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
    # Ignore first FLINK printed as test of configuration reloading
    elif line.startswith('START '):
      starting = 1
    elif line.startswith('FLINK '):
      if fresh == 1 and starting == 1:
        starting = 0
        link_width = float('nan')
      else:
        temp = line.split()
        ave_link = float(temp[3])
        print >> FLINK, "%g,%g" % (MDTU, ave_link)
        if len(temp) > 4:
          link_width = float(temp[4])   # To be printed with other widths
        else:
          link_width = float('nan')
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # plaquette, Polyakov loop, bosonic action, energy and CG iterations
    # Print both bosonic action and energy since lambda non-trivial
    # First four need Nc and DIMF extracted above
    elif line.startswith('GMES '):
      temp = line.split()
      print >> PLAQ, "%g,%g" % (MDTU, float(temp[4]) / Nc)
      print >> CG_ITERS, "%g,%g" % (traj, float(temp[3]))
      sB = float(temp[5]) / (1.5 * DIMF)
      print >> SB, "%g,%g" % (MDTU, sB)
      print >> ENERGY, "%g,%g" % (MDTU, 3.0 * (1.0 - sB) / lam)

      poly_r = float(temp[1]) / Nc
      poly_i = float(temp[2]) / Nc
      print >> POLY, "%g,%g" % (poly_r, poly_i)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD, "%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i)

      # Hack to account for failed polar decomposition
      # if first traj(s) after fresh start are rejected
      if fresh < 0 and poly_r == 1 and poly_i == 0:
        print >> POLY_MOD_POLAR, "%g,1,1,0" % MDTU
        print >> LINES_MOD_POLAR, "%g,1,1,1,1" % MDTU
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Wilson line in spatial direction
    elif line.startswith('LINES '):
      NEED_LINES = -1
      temp = line.split()
      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      print >> LINES, "%g,%g" % (x_r, x_i)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      print >> LINES_MOD, "%g,%g" % (MDTU, x_mod)

    # Unitarized Polyakov loop and spatial Wilson line
    elif line.startswith('LINES_POLAR '):
      temp = line.split()
      poly_r = float(temp[3]) / Nc
      poly_i = float(temp[4]) / Nc
      print >> POLY_POLAR, "%g,%g" % (poly_r, poly_i)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD_POLAR, "%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i)

      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      print >> LINES_POLAR, "%g,%g" % (x_r, x_i)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      print >> LINES_MOD_POLAR, "%g,%g" % (MDTU, x_mod)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Plaquette determinant and widths
    # Originally only in measurements, now can be in main output
    elif line.startswith('DET '):
      NEED_DET = -1
      temp = line.split()
      det_r = float(temp[1])
      det_i = float(temp[2])
      # !!! Site-by-site |1-det|^2 not measured on all ensembles
      # If it's missing, for now use volume average instead
      if len(temp) == 6:
        det = float(temp[5])
      else:
        det = (1.0 - det_r)**2 + det_i**2
      print >> DET, "%g,%g,%g,%g" % (MDTU, det, 1.0 - det_r, det_i)

    elif line.startswith('WIDTHS '):
      NEED_WIDTHS = -1
      temp = line.split()
      plaq_width = float(temp[1])
      re_width = float(temp[2])
      im_width = float(temp[3])
      print >> WIDTHS, "%g,%g,%g,%g,%g" \
                       % (MDTU, plaq_width, re_width, im_width, link_width)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Scalar eigenvalues
    # Only look at largest and smallest (most negative) to keep this manageable
    elif line.startswith('POLAR_EIG '):
      NEED_SCALAR_EIGS = -1
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
  # These are always paired (checked by check_eig_pairs.py)
  # Focus on first six pairs, 0, 2, 4, 6, 8 and 10
  infile = 'Out/eig.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
  else:
    check = -1    # Check whether file completed successfully
    temp = float('nan')
    eig = [temp, temp, temp, temp, temp, temp]
    big = [temp, temp, temp, temp, temp, temp]

    # Go
    for line in open(infile):
      if line.startswith('Time stamp '):
        stamp = line.rstrip()
        if stamp != oldstamp:
          print infile, "time stamp doesn't match final", oldstamp
          print >> ERRFILE, infile, "time stamp doesn't match final", oldstamp

      elif line.startswith('EIGENVALUE '):
        temp = line.split()
        index = int(temp[1])
        dat = float(temp[2])
        if index < 11 and index % 2 == 0:
          eig[index / 2] = dat
        if index == 0 and dat < min_eig:        # Check spectral range
          print infile, "exceeds RHMC spectral range:",
          print "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)
          print >> ERRFILE, infile, "exceeds RHMC spectral range:",
          print >> ERRFILE, "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)

      elif line.startswith('BIGEIGVAL  '):
        temp = line.split()
        index = int(temp[1])
        dat = float(temp[2])
        if index < 11 and index % 2 == 0:
          big[index / 2] = dat
        if index == 0 and dat > max_eig:        # Check spectral range
          print infile, "exceeds RHMC spectral range:",
          print "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)
          print >> ERRFILE, infile, "exceeds RHMC spectral range:",
          print >> ERRFILE, "%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig)

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

    print >> EIG, "%g,%g,%g,%g,%g,%g,%g" \
                  % (MDTU, eig[0], eig[1], eig[2], eig[3], eig[4], eig[5])
    print >> BIGEIG, "%g,%g,%g,%g,%g,%g,%g" \
                     % (MDTU, big[0], big[1], big[2], big[3], big[4], big[5])

    # Monitor (log of) condition number
    cond_num = math.log(big[0] / eig[0])
    print >> COND_NUM, "%d,%g" % (traj, cond_num)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "meas" file, if it is present
  # Main targens are Tr[X^2] and fermion-bilinear Ward identity
  infile = 'Out/meas.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
    print >> ERRFILE, "Problem opening", infile
    continue

  # We have a file, so let's cycle over its lines
  check = -1
  for line in open(infile):
    if line.startswith('Time stamp '):
      stamp = line.rstrip()
      if stamp != oldstamp:
        print infile, "time stamp doesn't match final", oldstamp
        print >> ERRFILE, infile, "time stamp doesn't match final", oldstamp

    # Need to worry about C2 in gauge action for Ward identity
    elif line.startswith('C2='):
      temp = line.split('=')
      C2 = float(temp[1])

    elif line.startswith('FLINK '): # Will be printed with other widths
      temp = line.split()
      if len(temp) > 7:
        link_width = float(temp[7])
      else:
        link_width = float('nan')

    # Scalar square Tr[X^2]
    elif line.startswith('TR_XSQ '):
      temp = line.split()
      tr_xsq = float(temp[1])
      print >> SCALARSQ, "%g,%g" % (MDTU, tr_xsq )

    # ----------------------------------------------------------
    # Fermion bilinear Ward identity
    elif 'CONGRAD' in line:
      CG = -1
    elif line.startswith('SUSY '):
      temp = line.split()
      trace = float(temp[1])
      gauge = float(temp[3])
      a = C2 * gauge
      b = trace
      susy = (a - b) / math.sqrt(a * a + b * b)

      # The imaginary part of the bilinear should vanish on average,
      # but large fluctuations may signal pathology
      zero = float(temp[2])
      print >> BILIN, "%g,%g,%g" % (MDTU, susy, zero)
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Wilson lines in the spatial direction
    elif NEED_LINES > 0 and line.startswith('LINES '):
      temp = line.split()
      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      print >> LINES, "%g,%g" % (x_r, x_i)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      print >> LINES_MOD, "%g,%g" % (MDTU, x_mod)

    # Unitarized Polyakov loop and Wilson lines in other directions
    elif NEED_LINES > 0 and line.startswith('LINES_POLAR '):
      temp = line.split()
      poly_r = float(temp[3]) / Nc
      poly_i = float(temp[4]) / Nc
      print >> POLY_POLAR, "%g,%g" % (poly_r, poly_i)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD_POLAR, "%g,%g,%g,%g" \
                               % (MDTU, poly_mod, poly_r, poly_i)

      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      print >> LINES_POLAR, "%g,%g" % (x_r, x_i)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      print >> LINES_MOD_POLAR, "%g,%g" % (MDTU, x_mod)
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Plaquette determinant and widths
    elif NEED_DET > 0 and line.startswith('DET '):
      temp = line.split()
      det_r = float(temp[1])
      det_i = float(temp[2])
      # !!! Site-by-site |1-det|^2 not measured on all ensembles
      # If it's missing, for now use volume average instead
      if len(temp) == 6:
        det = float(temp[5])
      else:
        det = (1.0 - det_r)**2 + det_i**2
      print >> DET, "%g,%g,%g,%g" % (MDTU, det, 1.0 - det_r, det_i)

    elif NEED_WIDTHS > 0 and line.startswith('WIDTHS '):
      temp = line.split()
      plaq_width = float(temp[1])
      re_width = float(temp[2])
      im_width = float(temp[3])
      print >> WIDTHS, "%g,%g,%g,%g,%g" \
                       % (MDTU, plaq_width, re_width, im_width, link_width)
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Scalar eigenvalues
    # Only look at largest and smallest (most negative) to keep this manageable
    elif NEED_SCALAR_EIGS > 0 and line.startswith('POLAR_EIG '):
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
    # ---------------------------------------------------------

  elif line.startswith('RUNNING COMPLETED'):
    if check == 1:    # Check that we have one measurement per file
      print infile, "reports two measurements"
      print >> ERRFILE, infile, "reports two measurements"
    check = 1
  if CG == -1:
    print infile, "encountered CG non-convergence"
    print >> ERRFILE, infile, "encountered CG non-convergence"
    CG = 1
  if check == -1:
    print infile, "did not complete"
    print >> ERRFILE, infile, "did not complete"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ERRFILE.close()
MISSINGFILES.close()
PLAQ.close()
SB.close()
SF.close()
POLY.close()
POLY_POLAR.close()
POLY_MOD.close()
POLY_MOD_POLAR.close()
LINES.close()
LINES_POLAR.close()
LINES_MOD.close()
LINES_MOD_POLAR.close()
FLINK.close()
DET.close()
BILIN.close()
WIDTHS.close()
SCALAR_EIG_AVE.close()
SCALAR_EIG.close()
SCALAR_EIG_WIDTHS.close()
EIG.close()
BIGEIG.close()
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
ENERGY.close()
SCALARSQ.close()
# ------------------------------------------------------------------
