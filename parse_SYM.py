#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse N=4 SYM output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting
# Normalize the Polyakov loop data by Nc and the bosonic action by 9Nc^2 / 2

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print("ERROR: Out/ does not exist")
  sys.exit(1)

ERRFILE = open('ERRORS', 'w')
MISSINGFILES = open('MISSING', 'w')

# Physical observables
PLAQ = open('data/plaq.csv', 'w')
print("MDTU,plaq_ss,plaq_st", file=PLAQ)
SB = open('data/SB.csv', 'w')
print("MDTU,S_B", file=SB)
SF = open('data/SF.csv', 'w')
print("MDTU,S_F", file=SF)
POLY = open('data/poly.csv', 'w')
print("ReTr(L),ImTr(L)", file=POLY)
POLY_POLAR = open('data/poly_polar.csv', 'w')
print("ReTr(L),ImTr(L)", file=POLY_POLAR)
POLY_MOD = open('data/poly_mod.csv', 'w')
print("MDTU,|Tr(L)|,ReTr(L),ImTr(L)", file=POLY_MOD)
POLY_MOD_POLAR = open('data/poly_mod_polar.csv', 'w')
print("MDTU,|Tr(L)|,ReTr(L),ImTr(L)", file=POLY_MOD_POLAR)
LINES = open('data/lines.csv', 'w')
print("ReTr(L*),ImTr(Lx),ImTr(Ly),ImTr(Lz),ImTr(L5)", file=LINES)
LINES_POLAR = open('data/lines_polar.csv', 'w')
print("ReTr(L*),ImTr(Lx),ImTr(Ly),ImTr(Lz),ImTr(L5)", file=LINES_POLAR)
LINES_MOD = open('data/lines_mod.csv', 'w')
print("MDTU,|Tr(Lx)|,|Tr(Ly)|,|Tr(Lz)|,|Tr(L5)|", file=LINES_MOD)
LINES_MOD_POLAR = open('data/lines_mod_polar.csv', 'w')
print("MDTU,|Tr(Lx)|,|Tr(Ly)|,|Tr(Lz)|,|Tr(L5)|", file=LINES_MOD_POLAR)
FLINK = open('data/Flink.csv', 'w')
print("MDTU,link", file=FLINK)
DET = open('data/det.csv', 'w')
print("MDTU,|det-1|^2,1-Re(det),Im(det)", file=DET)
WIDTHS = open('data/widths.csv', 'w')
print("MDTU,plaq,Re(det),Im(det),link", file=WIDTHS)
SCALAR_EIG_AVE = open('data/scalar_eig_ave.csv', 'w')
print("MDTU,min_ave,max_ave", file=SCALAR_EIG_AVE)
SCALAR_EIG = open('data/scalar_eig.csv', 'w')
print("MDTU,min_min,min_max,max_min,max_max5", file=SCALAR_EIG)
SCALAR_EIG_WIDTHS = open('data/scalar_eig_widths.csv', 'w')
print("MDTU,min_width,max_width", file=SCALAR_EIG_WIDTHS)
BILIN = open('data/bilin.csv', 'w')
print("MDTU,susyTrans,Im(bilin)", file=BILIN)
MONO = open('data/mono.csv', 'w')
print("MDTU,rho_M", file=MONO)
EIG = open('data/eig.csv', 'w')
print("MDTU,0,2,4,6,8,10", file=EIG)
BIGEIG = open('data/bigeig.csv', 'w')
print("MDTU,0,2,4,6,8,10", file=BIGEIG)

# Only create Wilson line eigenvalue data files if Out/WLeig* files present
do_WLeig = len(glob.glob('Out/WLeig.*'))
if do_WLeig > 0:
  WL_EIG = open('data/WLeig.csv', 'w')
  print("MDTU,eig0,eig1,eig2,eig3,eig4,eig5,eig6,eig7", file=WL_EIG)

# Evolution observables
ACCP = open('data/accP.csv', 'w')
print("t,accP", file=ACCP)
EXP_DS = open('data/exp_dS.csv', 'w')
print("t,e^(-dS)", file=EXP_DS)
DELTAS = open('data/deltaS.csv', 'w')
print("t,deltaS", file=DELTAS)
ABS_DS = open('data/abs_dS.csv', 'w')
print("t,|deltaS|", file=ABS_DS)
FORCE = open('data/force.csv', 'w')
print("t,G,F1,F2,F3", file=FORCE)
CG_ITERS = open('data/cg_iters.csv', 'w')
print("t,cg_iters", file=CG_ITERS)
WALLTIME = open('data/walltime.csv', 'w')
print("t,walltime", file=WALLTIME)
WALLTU = open('data/wallTU.csv', 'w')
print("t,cost", file=WALLTU)
COND_NUM = open('data/cond_num.csv', 'w')
print("t,cond_num", file=COND_NUM)

# Run parameters
NSTEP = open('data/Nstep.csv', 'w')
print("t,N_f,N_g", file=NSTEP)
STEPSIZE = open('data/stepsize.csv', 'w')
print("t,eps_f,eps_g", file=STEPSIZE)
NORDER = open('data/Norder.csv', 'w')
print("t,Norder", file=NORDER)
TLENGTH = open('data/tlength.csv', 'w')
print("t,L", file=TLENGTH)
KEY = open('data/key.csv', 'w')
print("t,file", file=KEY)
TU = open('data/TU.csv', 'w')
print("t,MDTU", file=TU)

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
    print("No Out/out.* files found")
    print("No Out/out.* files found", file=ERRFILE)
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
    print("Problem opening", infile)
    print("Problem opening", infile, file=ERRFILE)
    continue    # Skip to next file
#  print infile  # Monitor running status

  # If not starting from first file in this ensemble,
  # or if we seem to have skipped a file,
  # guess approximate starting trajectory
  traj_per_file = -1
  vol = -1
  for line in open(infile):
    if line.startswith('PLACEHOLDER'):
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      walltime = -2

    # Extract Nc for bosonic action and Polyakov loop normalizations
    # Convert it to DIMF to handle SU(N) runs
    # Should no longer need to handle pre-2014 formatting
    elif line.startswith('N=4 SYM, '):
      temp = line.split(',')
      Nc = float(((temp[1]).split())[2])
      DIMF = Nc**2
      temp = os.getcwd()
      if 'slnc' in temp:
        DIMF -= 1

    # Extract volume for normalizations
    elif line.startswith('nx '):
      vol = float((line.split())[1])
    elif line.startswith('ny ') or line.startswith('nz '):
      vol *= float((line.split())[1])
    elif line.startswith('nt '):
      vol *= float((line.split())[1])

    # Make sure no overflow errors parsing seed
    elif line.startswith('iseed '):
      if int((line.split())[1]) == -1:
        print(infile, "has suspicious seed -1")
        print(infile, "has suspicious seed -1", file=ERRFILE)

    elif line.startswith('trajecs '):
      traj_per_file = int((line.split())[1])
      endtraj = traj + traj_per_file
      break       # Don't go through whole file yet

  if traj_per_file < 0:
    print(infile, "never defines number of trajectories")
    print(infile, "never defines number of trajectories", file=ERRFILE)
    continue    # Skip to next file
  elif vol < 0:
    print(infile, "never defines lattice volume")
    print(infile, "never defines lattice volume", file=ERRFILE)
    continue    # Skip to next file

  if (traj == 0 and int(load) > 0) or (int(load) != oldcfg):
    print(infile, "misses MDTU before", load)
    traj = int(load)
    endtraj = traj + traj_per_file
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Cycle through lines in the file
  goodtogo = True
  for line in open(infile):
    # Check for premature termination (e.g., layout problem)
    if line.startswith('termination'):
      print(infile, "reports premature termination")
      print(infile, "reports premature termination", file=ERRFILE)
      goodtogo = False
      break

    # Check for unitarity problem
    # Report for info, but don't skip if found
    if line.startswith('Unitarity problem'):
      print(infile, "reports unitarity problem")
      print(infile, "reports unitarity problem", file=ERRFILE)
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
        print("%d,%d" % (endtraj, Norder), file=NORDER)
      else:         # Original 15-pole format didn't state spectral range
        min_eig = 1.0e-7
        max_eig = 1000.0
        print("%d,%d" % (endtraj, 15), file=NORDER)

    # Extract constant run parameters
    elif line.startswith('traj_length '):
      tlength = float((line.split())[1])
      print("%d,%g" % (endtraj, tlength), file=TLENGTH)
    elif line.startswith('nstep '):
      Nstep = int((line.split())[1])
      stepsize = tlength / float(Nstep)
    elif line.startswith('nstep_gauge '):
      Nstep_gauge = int((line.split())[1])
      stepsize_gauge = stepsize / float(2.0 * Nstep_gauge)
      Nstep_gauge *= 2 * Nstep
      print("%d,%d,%d" % (endtraj, Nstep, Nstep_gauge), file=NSTEP)
      print("%d,%g,%g" % (endtraj, stepsize, stepsize_gauge), file=STEPSIZE)

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
          print(infile, "time stamp doesn't match final", oldstamp)
          print(infile, "time stamp doesn't match final", \
                oldstamp, file=ERRFILE)
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
        fermAct[0] = float((line.split())[8]) / (16.0 * vol * DIMF)
      elif fermAct[1] < 0:
        fermAct[1] = float((line.split())[8]) / (16.0 * vol * DIMF)
      else:
        print(infile, "lists too many action computations")
        print(infile, "lists too many action computations", file=ERRFILE)

    elif ' delta S = ' in line:
      traj += 1
      temp = MDTU + tlength
      MDTU = round(temp, 3)   # Round off digits in MDTU
      print("%d,%g" % (traj, MDTU), file=TU)
      print("%g,%s" % (traj, tag), file=KEY)

      dS = (line.split())[4]
      # Strip nasty comma from dS, if necessary
      if ',' in dS:
        dS = (dS.split(','))[0]
      print("%d,%g" % (traj, float(dS)), file=DELTAS)
      print("%d,%g" % (traj, math.exp(-1.0 * float(dS))), file=EXP_DS)

      # For RMS, don't have an easy way to average over many measurements
      # Instead just print out absolute value and consider its running average
      print("%d,%g" % (traj, abs(float(dS))), file=ABS_DS)

      # Acceptance is smeared out by running averages
      if line.startswith('ACCEPT'):
        print("%d,1" % traj, file=ACCP)
        print("%d,%g" % (MDTU, fermAct[1]), file=SF)   # New action
      else:
        print("%d,0" % traj, file=ACCP)
        print("%d,%g" % (MDTU, fermAct[0]), file=SF)   # Original action
      fermAct = [-1.0, -1.0]                        # Reset

    # Forces -- take maxima rather than average if possible
    elif line.startswith('MONITOR_FORCE_GAUGE '):
      force_g = float((line.split())[-1])
    elif line.startswith('MONITOR_FORCE_FERMION0 '):
      force_f = float((line.split())[-1])
      if Nroot == 1:
        print("%d,%g,%g,null,null" % (traj, force_g, force_f), file=FORCE)
    elif line.startswith('MONITOR_FORCE_FERMION1 '):
      force_f2 = float((line.split())[-1])
      if Nroot == 2:
        print("%d,%g,%g,%g,null" % (traj, force_g, force_f, force_f2), \
              file=FORCE)
    elif line.startswith('MONITOR_FORCE_FERMION2 '):
      force_f3 = float((line.split())[-1])
      if Nroot == 3:
        print("%d,%g,%g,%g,%g" \
              % (traj, force_g, force_f, force_f2, force_f3), file=FORCE)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Ignore first FLINK printed as test of configuration reloading
    elif line.startswith('START '):
      starting = 1
    elif line.startswith('FLINK '):
      if fresh < 0 and starting == 1:
        starting = 0
        link_width = float('nan')
      else:
        temp = line.split()
        ave_link = float(temp[6])
        print("%g,%g" % (MDTU, ave_link), file=FLINK)
        if len(temp) > 7:
          link_width = float(temp[7]) # To be printed with other widths
        else:
          link_width = float('nan')
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # plaquette, Polyakov loop, bosonic action... and CG iterations
    # Normalize first three using Nc extracted above
    elif line.startswith('GMES '):
      temp = line.split()
      ss_plaq = float(temp[4]) / Nc
      st_plaq = float(temp[5]) / Nc
      print("%g,%g,%g" % (MDTU, ss_plaq, st_plaq), file=PLAQ)
      print("%g,%g" % (traj, float(temp[3])), file=CG_ITERS)

      print("%g,%g" % (MDTU, float(temp[6]) / (4.5 * DIMF)), file=SB)

      poly_r = float(temp[1]) / Nc
      poly_i = float(temp[2]) / Nc
      print("%g,%g" % (poly_r, poly_i), file=POLY)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print("%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i), file=POLY_MOD)

      # Hack to account for failed polar decomposition
      # if first traj(s) after fresh start are rejected
      if fresh == 1 and poly_r == 1 and poly_i == 0:
        print("%g,1,1,0" % MDTU, file=POLY_MOD_POLAR)
        print("%g,1,1,1,1" % MDTU, file=LINES_MOD_POLAR)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Wilson lines in other directions
    elif line.startswith('LINES '):
      NEED_LINES = -1
      temp = line.split()
      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      y_r = float(temp[3]) / Nc
      y_i = float(temp[4]) / Nc
      z_r = float(temp[5]) / Nc
      z_i = float(temp[6]) / Nc
      f_r = float(temp[9]) / Nc
      f_i = float(temp[10]) / Nc
      print("%g,%g,null,null,null" % (x_r, x_i), file=LINES)
      print("%g,null,%g,null,null" % (y_r, y_i), file=LINES)
      print("%g,null,null,%g,null" % (z_r, z_i), file=LINES)
      print("%g,null,null,null,%g" % (f_r, f_i), file=LINES)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      y_mod = math.sqrt(y_r**2 + y_i**2)
      z_mod = math.sqrt(z_r**2 + z_i**2)
      f_mod = math.sqrt(f_r**2 + f_i**2)
      print("%g,%g,%g,%g,%g" % (MDTU, x_mod, y_mod, z_mod, f_mod), \
            file=LINES_MOD)

    # Unitarized Polyakov loop and Wilson lines in other directions
    elif line.startswith('LINES_POLAR '):
      temp = line.split()
      poly_r = float(temp[7]) / Nc
      poly_i = float(temp[8]) / Nc
      print("%g,%g" % (poly_r, poly_i), file=POLY_POLAR)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print("%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i), \
            file=POLY_MOD_POLAR)

      x_r = float(temp[1]) / Nc
      x_i = float(temp[2]) / Nc
      y_r = float(temp[3]) / Nc
      y_i = float(temp[4]) / Nc
      z_r = float(temp[5]) / Nc
      z_i = float(temp[6]) / Nc
      f_r = float(temp[9]) / Nc
      f_i = float(temp[10]) / Nc
      print("%g,%g,null,null,null" % (x_r, x_i), file=LINES_POLAR)
      print("%g,null,%g,null,null" % (y_r, y_i), file=LINES_POLAR)
      print("%g,null,null,%g,null" % (z_r, z_i), file=LINES_POLAR)
      print("%g,null,null,null,%g" % (f_r, f_i), file=LINES_POLAR)
      x_mod = math.sqrt(x_r**2 + x_i**2)
      y_mod = math.sqrt(y_r**2 + y_i**2)
      z_mod = math.sqrt(z_r**2 + z_i**2)
      f_mod = math.sqrt(f_r**2 + f_i**2)
      print("%g,%g,%g,%g,%g" \
            % (MDTU, x_mod, y_mod, z_mod, f_mod), file=LINES_MOD_POLAR)
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
      print("%g,%g,%g,%g" % (MDTU, det, 1.0 - det_r, det_i), file=DET)

    elif line.startswith('WIDTHS '):
      NEED_WIDTHS = -1
      temp = line.split()
      plaq_width = float(temp[1])
      re_width = float(temp[2])
      im_width = float(temp[3])
      print("%g,%g,%g,%g,%g" \
            % (MDTU, plaq_width, re_width, im_width, link_width), file=WIDTHS)
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
        print("%g%s" % (MDTU, scalar_eig_ave), file=SCALAR_EIG_AVE)
        print("%g%s" % (MDTU, scalar_eig_ext), file=SCALAR_EIG)
        print("%g%s" % (MDTU, scalar_eig_width), file=SCALAR_EIG_WIDTHS)
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
    print(infile, "encountered CG non-convergence")
    print(infile, "encountered CG non-convergence", file=ERRFILE)
    CG = 1
  if walltime == -1:
    print(infile, "didn't print final timing")
    print(infile, "didn't print final timing", file=ERRFILE)
  elif walltime == -2:
    # Placeholder file -- error has been addressed as well as possible,
    # but don't print nonsense wall clock time
    fermAct = [-1.0, -1.0]                        # Reset
    pass
  else:   # We are good to go
    ave_time = walltime / traj_per_file
    print("%d,%g" % (traj, ave_time), file=WALLTIME)
    TUtime = walltime * float(cpus) / (60.0 * tlength * traj_per_file)
    print("%d,%g" % (traj, TUtime), file=WALLTU)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "eig" file, if it is present
  # These are always paired (checked by check_eig_pairs.py)
  # Focus on first six pairs, 0, 2, 4, 6, 8 and 10
  infile = 'Out/eig.' + cfg
  if not os.path.isfile(infile):
    print(infile, file=MISSINGFILES)
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
          print(infile, "time stamp doesn't match final", oldstamp)
          print(infile, "time stamp doesn't match final", oldstamp, \
                file=ERRFILE)

      elif line.startswith('EIGENVALUE '):
        temp = line.split()
        index = int(temp[1])
        dat = float(temp[2])
        if index < 11 and index % 2 == 0:
          eig[index / 2] = dat
        if index == 0 and dat < min_eig:        # Check spectral range
          print(infile, "exceeds RHMC spectral range:", end=' ')
          print("%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig))
          print(infile, "exceeds RHMC spectral range:", end=' ', file=ERRFILE)
          print("%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig), \
                file=ERRFILE)

      elif line.startswith('BIGEIGVAL  '):
        temp = line.split()
        index = int(temp[1])
        dat = float(temp[2])
        if index < 11 and index % 2 == 0:
          big[index / 2] = dat
        if index == 0 and dat > max_eig:        # Check spectral range
          print(infile, "exceeds RHMC spectral range:", end=' ')
          print("%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig))
          print(infile, "exceeds RHMC spectral range:", end=' ', file=ERRFILE)
          print("%.4g not in [%.4g, %.4g]" % (dat, min_eig, max_eig), \
                file=ERRFILE)

      elif 'WARNING' in line:
        print(infile, "saturated eigenvalue iterations")
        print(infile, "saturated eigenvalue iterations", file=ERRFILE)
      elif line.startswith('RUNNING COMPLETED'):
        if check == 1:    # Check that we have one measurement per file
          print(infile, "reports two measurements")
          print(infile, "reports two measurements", file=ERRFILE)
        check = 1
    if check == -1:
      print(infile, "did not complete")
      print(infile, "did not complete", file=ERRFILE)

    print("%g,%g,%g,%g,%g,%g,%g" \
          % (MDTU, eig[0], eig[1], eig[2], eig[3], eig[4], eig[5]), file=EIG)
    print("%g,%g,%g,%g,%g,%g,%g" \
          % (MDTU, big[0], big[1], big[2], big[3], big[4], big[5]), \
          file=BIGEIG)

    # Monitor (log of) condition number
    cond_num = math.log(big[0] / eig[0])
    print("%d,%g" % (traj, cond_num), file=COND_NUM)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "corr" file, if it is present
  # For now just care about fermion bilinears and related quantities
  infile = 'Out/corr.' + cfg
  if not os.path.isfile(infile):
    print(infile, file=MISSINGFILES)
  else:
    # Need to worry about C2 in gauge action
    # which wasn't always printed in output (though it is now)
    # For now, extract it from the path
    C2 = 1.0
    temp = os.getcwd()
    if '-c' in temp:
      temp2 = temp.split('-c')
      C2 = float(((temp2[1]).split('/'))[0])

    # We have a file, so let's cycle over its lines
    check = -1
    for line in open(infile):
      if line.startswith('Time stamp '):
        stamp = line.rstrip()
        if stamp != oldstamp:
          print(infile, "time stamp doesn't match final", oldstamp)
          print(infile, "time stamp doesn't match final", oldstamp, \
                file=ERRFILE)

      elif line.startswith('FLINK '): # Will be printed with other widths
        temp = line.split()
        if len(temp) > 7:
          link_width = float(temp[7])
        else:
          link_width = float('nan')

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
        print("%g,%g,%g" % (MDTU, susy, zero), file=BILIN)
      # ----------------------------------------------------------

      # ----------------------------------------------------------
      # Monopole world line density
      elif 'WARNING: total_mono mismatch' in line:
        print(infile, "has total_mono mismatch")
        print(infile, "has total_mono mismatch", file=ERRFILE)
      elif line.startswith('MONOPOLE '):
        mono = float((line.split())[10])
        print("%g,%g" % (MDTU, mono / (4.0 * vol)), file=MONO)
      # ----------------------------------------------------------

      # ----------------------------------------------------------
      # Wilson lines in other directions
      elif NEED_LINES > 0 and line.startswith('LINES '):
        temp = line.split()
        x_r = float(temp[1]) / Nc
        x_i = float(temp[2]) / Nc
        y_r = float(temp[3]) / Nc
        y_i = float(temp[4]) / Nc
        z_r = float(temp[5]) / Nc
        z_i = float(temp[6]) / Nc
        f_r = float(temp[9]) / Nc
        f_i = float(temp[10]) / Nc
        print("%g,%g,null,null,null" % (x_r, x_i), file=LINES)
        print("%g,null,%g,null,null" % (y_r, y_i), file=LINES)
        print("%g,null,null,%g,null" % (z_r, z_i), file=LINES)
        print("%g,null,null,null,%g" % (f_r, f_i), file=LINES)
        x_mod = math.sqrt(x_r**2 + x_i**2)
        y_mod = math.sqrt(y_r**2 + y_i**2)
        z_mod = math.sqrt(z_r**2 + z_i**2)
        f_mod = math.sqrt(f_r**2 + f_i**2)
        print("%g,%g,%g,%g,%g" \
              % (MDTU, x_mod, y_mod, z_mod, f_mod), file=LINES_MOD)

      # Unitarized Polyakov loop and Wilson lines in other directions
      elif NEED_LINES > 0 and line.startswith('LINES_POLAR '):
        temp = line.split()
        poly_r = float(temp[7]) / Nc
        poly_i = float(temp[8]) / Nc
        print("%g,%g" % (poly_r, poly_i), file=POLY_POLAR)
        poly_mod = math.sqrt(poly_r**2 + poly_i**2)
        print("%g,%g,%g,%g" \
              % (MDTU, poly_mod, poly_r, poly_i), file=POLY_MOD_POLAR)

        x_r = float(temp[1]) / Nc
        x_i = float(temp[2]) / Nc
        y_r = float(temp[3]) / Nc
        y_i = float(temp[4]) / Nc
        z_r = float(temp[5]) / Nc
        z_i = float(temp[6]) / Nc
        f_r = float(temp[9]) / Nc
        f_i = float(temp[10]) / Nc
        print("%g,%g,null,null,null" % (x_r, x_i), file=LINES_POLAR)
        print("%g,null,%g,null,null" % (y_r, y_i), file=LINES_POLAR)
        print("%g,null,null,%g,null" % (z_r, z_i), file=LINES_POLAR)
        print("%g,null,null,null,%g" % (f_r, f_i), file=LINES_POLAR)
        x_mod = math.sqrt(x_r**2 + x_i**2)
        y_mod = math.sqrt(y_r**2 + y_i**2)
        z_mod = math.sqrt(z_r**2 + z_i**2)
        f_mod = math.sqrt(f_r**2 + f_i**2)
        print("%g,%g,%g,%g,%g" \
              % (MDTU, x_mod, y_mod, z_mod, f_mod), file=LINES_MOD_POLAR)
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
        print("%g,%g,%g,%g" % (MDTU, det, 1.0 - det_r, det_i), file=DET)

      elif NEED_WIDTHS > 0 and line.startswith('WIDTHS '):
        temp = line.split()
        plaq_width = float(temp[1])
        re_width = float(temp[2])
        im_width = float(temp[3])
        print("%g,%g,%g,%g,%g" \
              % (MDTU, plaq_width, re_width, im_width, link_width), \
              file=WIDTHS)
      # ----------------------------------------------------------

      # ----------------------------------------------------------
      # Scalar eigenvalues
      # Only look at largest and smallest (most negative)
      # to keep this manageable
      elif NEED_SCALAR_EIGS > 0 and line.startswith('POLAR_EIG '):
        temp = line.split()
        index = int(temp[1])
        if index == 0 or index == Nc - 1:
          scalar_eig_ave += ',' + str(temp[2])
          scalar_eig_ext += ',' + str(temp[4]) + ',' + str(temp[5])
          scalar_eig_width += ',' + str(temp[3])
        if index == Nc - 1:
          print("%g%s" % (MDTU, scalar_eig_ave), file=SCALAR_EIG_AVE)
          print("%g%s" % (MDTU, scalar_eig_ext), file=SCALAR_EIG)
          print("%g%s" % (MDTU, scalar_eig_width), file=SCALAR_EIG_WIDTHS)
          scalar_eig_ave = ''
          scalar_eig_ext = ''
          scalar_eig_width = ''
      # ---------------------------------------------------------



      elif line.startswith('RUNNING COMPLETED'):
        if check == 1:    # Check that we have one measurement per file
          print(infile, "reports two measurements")
          print(infile, "reports two measurements", file=ERRFILE)
        check = 1
    if CG == -1:
      print(infile, "encountered CG non-convergence")
      print(infile, "encountered CG non-convergence", file=ERRFILE)
      CG = 1
    if check == -1:
      print(infile, "did not complete")
      print(infile, "did not complete", file=ERRFILE)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "WLeig" file, if it is present
  # Only running on thermalized configurations, so don't report any missing
  infile = 'Out/WLeig.' + cfg
  if os.path.isfile(infile):
    # We have a file, so let's cycle over its lines
    check = -1
    for line in open(infile):
      if line.startswith('Time stamp '):
        stamp = line.rstrip()
        if stamp != oldstamp:
          print(infile, "time stamp doesn't match final", oldstamp)
          print(infile, "time stamp doesn't match final", oldstamp, \
                file=ERRFILE)

      # Format: LINES_POLAR_EIG x y z t dir {Nc x phase}
      # Check that all phases are within [-pi, pi),
      # accounting for rounding in output files
      # Can be commented out to speed up analysis
      elif line.startswith('LINES_POLAR_EIG '):
        WL_EIG.write("%g," % (MDTU))
        temp = line.split()
        for i in range(int(Nc - 1)):
          phase = float(temp[6 + i])
          if phase > 3.142 or phase < -3.142:
            print(infile, "phase %.4g exceeds [-pi, pi)" % (phase))
            print(infile, "phase %.4g exceeds [-pi, pi)" % (phase), \
                  file=ERRFILE)
          WL_EIG.write("%g," % (phase))
        index = int(5 + Nc)
        print("%g" % (float(temp[index])), file=WL_EIG)

      elif line.startswith('RUNNING COMPLETED'):
        if check == 1:    # Check that we have one measurement per file
          print(infile, "reports two measurements")
          print(infile, "reports two measurements", file=ERRFILE)
        check = 1
    if check == -1:
      print(infile, "did not complete")
      print(infile, "did not complete", file=ERRFILE)
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
MONO.close()
WIDTHS.close()
SCALAR_EIG_AVE.close()
SCALAR_EIG.close()
SCALAR_EIG_WIDTHS.close()
EIG.close()
BIGEIG.close()
if do_WLeig > 0:
 WL_EIG.close()

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
