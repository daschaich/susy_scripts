#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse N=(2, 2) SYM output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting
# Normalize the Polyakov loop data by Nc and the bosonic action by 3Nc^2 / 2

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
POLY = open('data/poly.csv', 'w')
print >> POLY, "ReTr(L),ImTr(L)"
POLY_MOD = open('data/poly_mod.csv', 'w')
print >> POLY_MOD, "MDTU,|Tr(L)|,ReTr(L),ImTr(L)"
FLINK = open('data/Flink.csv', 'w')
print >> FLINK, "MDTU,link"
DET = open('data/det.csv', 'w')
print >> DET, "MDTU,|det-1|,|Re(det)-1|"
EIG = open('data/eig.csv', 'w')
print >> EIG, "MDTU,0,2,4,6,8,10"
BILIN = open('data/bilin.csv', 'w')
print >> BILIN, "MDTU,susyTrans"
MONO = open('data/mono.csv', 'w')
print >> MONO , "MDTU,rho_M"

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
print >> FORCE, "t,G,F"
CG_ITERS = open('data/cg_iters.csv', 'w')
print >> CG_ITERS, "t,cg_iters"
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
oldcfg = 0
oldstamp = "start"
traj = 0;
MDTU = 0;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through files out.$load-$save from list.txt
for temp_tag in open('list.txt'):
  tag = temp_tag.rstrip()
  load, cfg = tag.split('-')
  # Initialize running sums and set dummy walltime
  # If the walltime isn't overwritten, then the run died
  # or its output file is corrupted
  walltime = -1;
  stamp = "start";

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
  traj_per_file = -1;
  vol = -1;
  for line in open(infile):
    if line.startswith('PLACEHOLDER'):
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      walltime = -2;

    # Extract volume for monopole world line density
    elif line.startswith('nx '):
      vol = float((line.split())[1])
    elif line.startswith('nt '):
      vol *= float((line.split())[1])

    # Extract Nc for bosonic action and Polyakov loop normalizations
    elif line.startswith('N=(2, 2) SYM, '):
      temp1 = line.split(',')
      Nc = float(((temp1[2]).split())[2])
    elif line.startswith('trajecs'):
      traj_per_file = int((line.split())[1])
      endtraj = traj + traj_per_file
      break       # Don't go through whole file yet

  if traj_per_file < 0:
    print infile, "never defines number of trajectories"
    print >> ERRFILE, infile, "never defines number of trajectories"
    continue    # Skip to next file
  elif vol < 0:
    print infile, "never defines lattice volume"
    print >> ERRFILE, infile, "never defines lattice volume"
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
    # Acceptance comes before measurements
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

      # Will be smeared out by running averages
      if line.startswith('ACCEPT'):
        print >> ACCP, "%d,1" % traj
      else:
        print >> ACCP, "%d,0" % traj

    # Forces -- take maxima rather than average if possible
    elif line.startswith('MONITOR_FORCE_GAUGE '):
      force_g = float((line.split())[-1])
    elif line.startswith('MONITOR_FORCE_FERMION0 '):
      force_f = float((line.split())[-1])
      print >> FORCE, "%d,%g,%g" % (traj, force_g, force_f)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    elif line.startswith('START '):
      starting = 1
    elif line.startswith('FLINK '):
      if starting == 1:
        starting = 0
      else:
        ave_link = float((line.split())[3])
        print >> FLINK, "%g,%g" % (MDTU, ave_link)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Gauge measurements come next:
    # plaquette, Polyakov loop, bosonic action... and CG iterations
    elif line.startswith('GMES '):
      temp = line.split()
      print >> PLAQ, "%g,%g" % (MDTU, float(temp[4]))
      print >> CG_ITERS, "%g,%g" % (traj, float(temp[3]))

      # Normalize bosonic action and Polyakov loop using Nc extracted above
      print >> SB, "%g,%g" % (MDTU, float(temp[5]) / (1.5 * Nc**2))

      poly_r = float(temp[1]) / Nc
      poly_i = float(temp[2]) / Nc
      print >> POLY, "%g,%g" % (poly_r, poly_i)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD, "%g,%g,%g,%g" % (MDTU, poly_mod, poly_r, poly_i)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Finally, the determinant is only measured every once in a while
    elif line.startswith('DET '):
      temp = line.split()
      det_r = abs(float(temp[1]) - 1.0)
      det_i = float(temp[2])
      det = math.sqrt(det_r**2 + det_i**2)
      print >> DET, "%g,%g,%g" % (MDTU, det, det_r)

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
    pass
  else:   # We are good to go
    ave_time = walltime / traj_per_file
    print >> WALLTIME, "%d,%g" % (traj, ave_time)
    TUtime = walltime * float(cpus) / (60.0 * tlength * traj_per_file)
    print >> WALLTU, "%d,%g" % (traj, TUtime)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "eig" file, if it is present
  # These are always paired (checked by check_eig_pairs.py),
  # Focus on first six pairs, 0, 2, 4, 6, 8 and 10
  infile = 'Out/eig.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
  else:
    check = -1    # Check whether file completed successfully
    temp = float('nan')
    eig = [temp, temp, temp, temp, temp, temp]

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
        if index < 11 and index % 2 == 0:
          eig[index / 2] = float(temp[2])
      elif 'WARNING' in line:
        print infile, "saturated eigenvalue iterations"
        print >> ERRFILE, infile, "saturated eigenvalue iterations"
      elif line.startswith('RUNNING COMPLETED'):
        check = 1
    if check == -1:
      print infile, "did not complete"
      print >> ERRFILE, infile, "did not complete"

    print >> EIG, "%g,%g,%g,%g,%g,%g,%g" \
                  % (MDTU, eig[0], eig[1], eig[2], eig[3], eig[4], eig[5])
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "corr" file, if it is present
  # For now just care about fermion bilinears and related quantities
  infile = 'Out/corr.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
  else:
    # Need to worry about C2 in gauge action
    # which wasn't always printed in output (though it is now)
    # For now, extract it from the path
    C2 = 1.0
    if '-c' in os.getcwd():
      temp1 = os.getcwd()
      temp2 = temp1.split('-c')
      C2 = float(((temp1[1]).split('/'))[0])

    # We have a file, so let's cycle over its lines
    check = -1
    for line in open(infile):
      if line.startswith('Time stamp '):
        stamp = line.rstrip()
        if stamp != oldstamp:
          print infile, "time stamp doesn't match final", oldstamp
          print >> ERRFILE, infile, "time stamp doesn't match final", oldstamp

      elif line.startswith('SUSY '):
        temp = line.split()
        trace = float(temp[1])
        gauge = float(temp[3])
        susy = (C2 * gauge - trace) / (C2 * gauge + trace)
        print >> BILIN, "%g,%g" % (MDTU, susy)
      elif 'WARNING' in line:
        print infile, "has total_mono mismatch"
        print >> ERRFILE, infile, "has total_mono mismatch"
      elif line.startswith('MONOPOLE'):
        mono = float((line.split())[6])
        print >> MONO, "%g,%g" % (MDTU, mono / (2.0 * vol))
      elif line.startswith('RUNNING COMPLETED'):
        check = 1
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
POLY.close()
POLY_MOD.close()
FLINK.close()
DET.close()
EIG.close()
BILIN.close()
MONO.close()
ACCP.close()
EXP_DS.close()
DELTAS.close()
ABS_DS.close()
FORCE.close()
CG_ITERS.close()
WALLTIME.close()
WALLTU.close()
NSTEP.close()
STEPSIZE.close()
TLENGTH.close()
KEY.close()
TU.close()
# ------------------------------------------------------------------
