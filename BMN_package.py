#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package BMN data and results/attributes into HDF5 file

# Cycle over ensembles and write to ~/SYM/BMN/BMN_data.h5
# Total of 292 ensembles: {31, 121, 104, 36} for Nc={4, 8, 12, 16}, resp.
# Group paths will specify Nc, Nt, g=lambda/mu^3 and f=T/mu
# Attributes for each ensemble:
#   mu_lat, lambda_lat, number of trajectories, thermalization cut, block size,
#   number of blocks, acceptance rate, and three autocorrelation times:
#     for |PL|, lowest fermion operator eigenvalue, and Tr[X_9^2]
# Datasets for each ensemble, most with ave and err as attributes:
#   bosonic action, ratio of SO(6) and SO(3) Tr[X^2] action contributions,
#   Myers term in the action, exp(-Delta S), extremal scalar eigenvalues,
#   real, imag and magnitude of Polyakov loop (PL)
#     (magnitude with susceptibility as attribute),
#   all scalar squares Tr[X^2], PL eigenvalue phases,
#   extremal fermion operator eigenvalues, condition number, spectral range,
#   list of configurations for which eigenvalues computed,
#   pfaffian phase and list of configs for which it was measured
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Helper function to map rational approximation order to its spectral range
def spect_range(N):
  if N == 5:
    return (0.1, 50)
  elif N == 6:
    return (0.02, 50)
  elif N == 7:
    return (0.01, 150)
  elif N == 8:
    return (0.001, 50)
  elif N == 9:
    return (1e-4, 45)
  elif N == 10:
    return (5e-4, 1000)
  elif N == 11:
    return (1e-5, 50)
  elif N == 12:
    return (5e-5, 2500)
  elif N == 13:
    return (5e-5, 5000)
  elif N == 14:
    return (1e-6, 1900)
  elif N == 15:
    return (1e-7, 1000)
  elif N == 16:
    return (1e-8, 500)
  elif N == 17:
    return (5e-8, 2500)
  elif N == 18:
    return (1e-9, 1000)
  elif N == 19:
    return (1e-8, 50000)
  else:
    print("ERROR: Unrecognized Norder %d" % Norder)
    sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Annoyingly, need to use absolute Barkla path
path = "/users/schaich/SYM/BMN/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)
f = h5py.File("/users/schaich/SYM/BMN/BMN_data.h5", 'w')

# Top-level groups for each Nc
for Nc in ['Nc4', 'Nc8', 'Nc12', 'Nc16']:
  f.create_group(Nc)

  # Second-level groups for each Nt
  if Nc == 'Nc4':
    allNt = ['8', '12', '16', '24', '32', '48', '64', '96', '128']
  elif Nc == 'Nc8' or Nc == 'Nc12':
    allNt = ['8', '16', '24']
  elif Nc == 'Nc16':
    allNt = ['16', '24']
  else:
    print("ERROR: Unrecognized Nc=%s" % Nc)
    sys.exit(1)

  for Nt in allNt:
    Nt_tag = 'Nt' + Nt
    top_dir = Nc + '_' + 'nt' + Nt + '/'
    this_Nc_Nt = Nc + '/' + Nt_tag
    f.create_group(this_Nc_Nt)

    # Third- and fourth-level groups for each (g, f=T/mu) ensemble
    os.chdir(path + top_dir)
    if top_dir == 'Nc8_nt8/':
      ens_list = glob.glob('g0.000*[0-9]') \
               + glob.glob('g0.001_f0.0[789]') \
               + glob.glob('g0.001_f0.1')
    elif top_dir == 'Nc12_nt8/':
      ens_list = glob.glob('g0.000*[0-9]')
    else:
      ens_list = glob.glob('g0.0*[0-9]')
    for ens in ens_list:
      os.chdir(path + top_dir + ens)
      toCheck = 'results/poly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished ensembles
        print("Skipping %s" % top_dir + ens)
        continue

      temp = ens.split('_')
      g = temp[0]
      T_ov_mu = temp[1]
      this_ens = this_Nc_Nt + '/' + g + '/' + T_ov_mu
      this_grp = f.create_group(this_ens)

      # Record mu_lat = 1 / (Nt f) and lambda_lat = g mu_lat^3
      Nt_float = float(Nt_tag.split('Nt')[1])
      g_float = float(g.split('g')[1])
      f_float = float(T_ov_mu.split('f')[1])
      mu_lat = 1.0 / (Nt_float * f_float)
      this_grp.attrs['mu_lat'] = mu_lat
      this_grp.attrs['lambda_lat'] = g_float * mu_lat**3

      # Record thermalization cut and block size
      check = -1
      therm = path + top_dir + 'therm.sh'
      for line in open(therm):
        if 'g0.0' in line:
          temp = line.split()
          if ens == temp[1]:
            this_grp.attrs['thermalization cut'] = int(temp[2])
            this_grp.attrs['block size'] = int(temp[3])
            check = int(temp[2])
            break
      if check < 0:
        print("ERROR: Thermalization cut not found for %s: " % this_str)
        sys.exit(1)

      # ------------------------------------------------------------
      # Record acceptance, auto-correlation and Nblocks for ensemble
      for line in open('results/accP.dat'):
        temp = line.split()
        Nblocks = int(temp[-1])
        this_grp.attrs['Nblocks'] = Nblocks
        this_grp.attrs['acceptance'] = float(temp[0])
      for line in open('results/poly_mod.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time (|PL|)'] = float(temp[0])
      for line in open('results/eig.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time (eig)'] = float(temp[0])
      for line in open('results/scalarsquares.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time (Tr[X^2])'] = float(temp[0])

      # Now set up datasets and results-attributes for this ensemble
      # First do bosonic action and set Ntraj
      SB_arr = []
      for line in open('data/SB.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        SB_arr.append(float(temp[1]))
      Ntraj = len(SB_arr)
      this_grp.attrs['Ntraj'] = Ntraj

      SB = np.array(SB_arr)
      dset = this_grp.create_dataset('SB', data=SB)

      # Record ensemble average as an attribute of the dataset
      for line in open('results/SB.dat'):
        temp = line.split()
        dset.attrs['ave'] = float(temp[0])
        dset.attrs['err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in plaq.suscept" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Simple observables measured every trajectory=MDTU
      # SO(6)/SO(3) ratio, Myers term, exp(-Delta S)
      # Only want one datum per line from each of these
      # Myers normalized while parsing
      # Skipping SF, cg_iters, wall_time
      for obs in ['ratio', 'Myers', 'exp_dS']:
        obsfile = 'data/' + obs + '.csv'
        dset = this_grp.create_dataset(obs, (Ntraj,), dtype='f')

        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            continue
          dset[traj] = float(temp[1])
          traj += 1

        if traj != Ntraj:
          print("ERROR: Ntraj mismatch in %s, " % this_str, end='')
          print("%d vs %d in %s" % (traj, Ntraj, obsfile))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        if os.path.isfile(resfile):           # No ave for poly_arg
          for line in open(resfile):
            if line.startswith('#'):
              continue
            temp = line.split()
            dset.attrs['ave'] = float(temp[0])
            dset.attrs['err'] = float(temp[1])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s, " % this_str, end='')
              print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
              sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Extremal scalar eigenvalues also measured every trajectory=MDTU
      # Want two data per line (average minimum and maximum eigenvalues)
      # These are not averaged
      dset = this_grp.create_dataset('scalar_eig_ave', (Ntraj,2,), dtype='f')
      dset.attrs['columns'] = ['min', 'max']
      traj = 0
      for line in open('data/scalar_eig_ave.csv'):
        temp = line.split(',')
        if line.startswith('M') or line.startswith('t'):
          continue
        dset[traj] = [float(temp[1]), float(temp[2])]
        traj += 1

      if traj != Ntraj:
        print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
        print("%d vs %d in scalar_eig_ave" % (traj, Ntraj))
        sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Polyakov loop also measured every trajectory=MDTU
      # Want three data per line from each of these (real, imag, magnitude)
      # Only magnitude is averaged
      obsfile = 'data/poly_mod.csv'
      name = 'Polyakov_loop'

      dset = this_grp.create_dataset(name, (Ntraj,3,), dtype='f')
      dset.attrs['columns'] = ['real', 'imag', 'magnitude']
      traj = 0
      for line in open(obsfile):
        temp = line.split(',')
        if line.startswith('M') or line.startswith('t'):
          continue
        dset[traj] = [float(temp[2]), float(temp[3]), float(temp[1])]
        traj += 1

      if traj != Ntraj:
        print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
        print("%d vs %d in %s" % (traj, Ntraj, obsfile))
        sys.exit(1)

      # Results (including susceptibility) as attributes, checking Nblocks
      for line in open('results/poly_mod.dat'):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['magnitude ave'] = float(temp[0])
        dset.attrs['magnitude err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
          print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
          sys.exit(1)
      for line in open('results/poly_mod.suscept'):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['magnitude suscept'] = float(temp[0])
        dset.attrs['magnitude suscept_err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
          print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Scalar squares Tr[X^2] also measured every trajectory=MDTU
      # Want nine data per line, all averaged
      # Here it is more convenient to collect the phases in a numpy array
      Nscalar = 9
      TrXSq = np.empty((Ntraj, Nscalar), dtype = np.float)

      TrXSq_tot = 0
      for line in open('data/scalarsquares.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        index = int(temp[0]) - 1
        for j in range(Nscalar):
          TrXSq[index][j] = float(temp[j + 1])
          TrXSq_tot += 1

      # Check to make sure no measurements are missing
      expect = Ntraj * Nscalar
      if TrXSq_tot != expect:
        print("ERROR: TrXSq mismatch in %s: " % this_ens, end='')
        print("%d vs %d expected" % (TrXSq_tot, expect))
        sys.exit(1)

      dset = this_grp.create_dataset('Scalar_squares (Tr[X^2])', data=TrXSq)

      # Results (including splitting) as attributes, checking Nblocks
      for line in open('results/scalarsquares.dat'):
        TrXSq_ave = np.empty((Nscalar), dtype = np.float)
        TrXSq_err = np.empty_like(TrXSq_ave)
        temp = line.split()
        for j in range(Nscalar):
          TrXSq_ave[j] = float(temp[2 * j])
          TrXSq_err[j] = float(temp[2 * j + 1])

        dset.attrs['ave'] = TrXSq_ave
        dset.attrs['err'] = TrXSq_err
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
          print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
          sys.exit(1)
      for line in open('results/splitting.dat'):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['SO(6) splitting ave'] = float(temp[0])
        dset.attrs['SO(6) splitting err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
          print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Polyakov loop eigenvalue phases
      # Nc phases measured every trajectory=MDTU, not averaged
      # Here it is more convenient to collect the phases in a numpy array
      Nc_int = int(Nc.split('Nc')[1])
      PLeig = np.empty((Ntraj, Nc_int), dtype = np.float)

      PLeig_tot = 0
      for line in open('data/PLeig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        index = int(temp[0]) - 1
        for j in range(Nc_int):
          PLeig[index][j] = float(temp[j + 1])
          PLeig_tot += 1

      # Check to make sure no measurements are missing
      expect = Ntraj * Nc_int
      if PLeig_tot != expect:
        print("ERROR: PLeig mismatch in %s: " % this_ens, end='')
        print("%d vs %d expected" % (PLeig_tot, expect))
        sys.exit(1)

      dset = this_grp.create_dataset('PLeig_phases', data=PLeig)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Fermion operator eigenvalues and RHMC spectral range
      # Format: min_eig, max_eig, log(cond_num), min_range, max_range
      # No averaging so attributes are only column labels
      # Save (and count) trajectories at which measurements were run
      meas_arr = []

      # First pass through to count number of measurements
      for line in open('data/eig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        meas_arr.append(int(temp[0]))
      Nmeas = len(meas_arr)
      eig = np.empty((Nmeas, 5), dtype = np.float)

      # Minimum eigenvalue from eig.csv
      eig_count = 0
      for line in open('data/eig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        eig[eig_count][0] = float(temp[1])
        eig_count += 1
      if eig_count != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in min_eig" % (meas, Nmeas))
        sys.exit(1)

      # Dataset for measurements
      meas = np.array(meas_arr)
      this_grp.create_dataset('Fermion_op_eig measurements', data=meas)

      # Maximum eigenvalue from bigeig.csv
      eig_count = 0
      for line in open('data/bigeig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        eig[eig_count][1] = float(temp[1])
        eig_count += 1
      if eig_count != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in max_eig" % (meas, Nmeas))
        sys.exit(1)

      # Logarithm of condition number
      eig_count = 0
      for line in open('data/cond_num.csv'):
        if line.startswith('t'):
          continue
        temp = line.split(',')
        eig[eig_count][2] = float(temp[1])
        eig_count += 1
      if eig_count != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in cond_num" % (meas, Nmeas))
        sys.exit(1)

      # Spectral range from order of rational approximation
      eig_count = 0
      for line in open('data/Norder.csv'):
        if line.startswith('t'):
          continue
        temp = line.split(',')
        srange = spect_range(int(temp[1]))
        eig[eig_count][3] = float(srange[0])
        eig[eig_count][4] = float(srange[1])
        eig_count += 1
      if eig_count != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in spectral range" % (meas, Nmeas))
        sys.exit(1)

      # Sanity checks on eigenvalue measurements
      for j in range(Nmeas):
        if eig[j][0] < eig[j][3] or eig[j][1] > eig[j][4]:
          print("ERROR: Spectral range problem in %s, " % this_ens, end='')
          print("[%.4g, %.4g] vs " % (eig[j][0], eig[j][1]), end='')
          print("[%.4g, %.4g] vs " % (eig[j][3], eig[j][4]), end='')
          sys.exit(1)

      dset = this_grp.create_dataset('Fermion_op_eig', data=eig)
      dset.attrs['columns'] = ['min_eig', 'max_eig', 'log(cond_num)', \
                               'min_range', 'max_range']
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Finally pfaffian phase measurements when present
      # (typically for small Nc and Nt)
      # Save (and count) configurations on which measurements were run
      pfile = 'results/pfaffian.dat'
      if os.path.isfile(pfile):        # Not present for most groups
        pfaff_arr = []
        conf_arr = []
        files = sorted(glob.glob('Out/phase.*'))
        Npfaff = len(files)
        for file in files:
          conf_arr.append(int(file.split('phase.')[1]))
          for line in open(file):
            # Format: PFAFF log_mag phase |cos(phase)| |sin(phase)|
            if line.startswith('PFAFF '):
              temp = line.split()
              pfaff_arr.append(float(temp[2]))

        conf = np.array(conf_arr)
        this_grp.create_dataset('Pfaffian measurements', data=conf)
        pfaff = np.array(pfaff_arr)
        dset = this_grp.create_dataset('Pfaffian phase', data=pfaff)
        
        # Check that all files have a measurement
        if not len(conf) == Npfaff:
          print("ERROR: Npfaff conf mismatch: %d vs %d in %s" \
                % (len(conf), Npfaff, top_dir + ens))
          sys.exit(1)
        if not len(pfaff) == Npfaff:
          print("ERROR: Npfaff phase mismatch: %d vs %d in %s" \
                % (len(pfaff), Npfaff, top_dir + ens))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        for line in open(pfile):
          if line.startswith('From'):
            temp = line.split()
            if not int(temp[1]) == Npfaff:
              print("ERROR: Npfaff ave mismatch: %d vs %d in %s" \
                    % (len(conf), Npfaff, top_dir + ens))
              sys.exit(1)
          elif line.startswith('Real: '):
            temp = line.split()
            dset.attrs['Real ave'] = float(temp[1])
            dset.attrs['Real err'] = float(temp[3].replace(";", ""))
            dset.attrs['Imag ave'] = float(temp[5])
            dset.attrs['Imag err'] = float(temp[7].replace(";", ""))
          elif line.startswith('Mag: '):
            temp = line.split()
            dset.attrs['Magnitude ave'] = float(temp[1])
            dset.attrs['Magnitude err'] = float(temp[3].replace(";", ""))
            dset.attrs['Phase ave'] = float(temp[5])
            dset.attrs['Phase err'] = float(temp[7].replace(";", ""))
# ------------------------------------------------------------------
