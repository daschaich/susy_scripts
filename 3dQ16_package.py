#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package 3d SYM data and results/attributes into HDF5 file

# Cycle over ensembles and write to ~/SYM/3d/3dSYM_data.h5
# 78 Nc=8 ensembles plus one Nc=4 12nt12 ensemble and one Nc=6 12nt12 ensemble
# Group paths will specify Nc, L=Nt,
#   rt = lambda_lat * Nt = 4 / (t * sqrt(3)),
#   and g = mu / lambda_lat
# Attributes for each ensemble:
#   t = 4 / (rt * sqrt(3)), lambda_lat, mu, kappa, number of trajectories,
#   thermalization cut, block size, number of blocks, acceptance rate,
#   autocorrelation times for |ML|, lowest eigenvalue,
#   and fermion bilinear Q-susy Ward identity violation
# Datasets for each ensemble, most with ave and err as attributes:
#   plaquette, bosonic action, Tr[U.Ubar]/N, exp(-Delta S),
#   real, imag and magnitude of Maldacena (ML) and Polyakov (PL) loops,
#   magnitudes of (complexified and unitarized) Wilson lines in spatial dirs,
#   fermion bilinear Q-susy Ward identify violation,
#   extremal eigenvalues, condition number, spectral range,
#   phases of Wilson line eigenvalues, lists of analyzed configurations
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
path = "/users/schaich/SYM/3d/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)
f = h5py.File("3dSYM_data.h5", 'w')

# Top-level groups for each Nc
for Nc in ['Nc4', 'Nc6', 'Nc8']:
  f.create_group(Nc)

  # Second- and third-level groups for each Nt and Nx=Ny=L
  # Currently we only have L=Nt, but let's be a bit flexible
  if Nc == 'Nc4' or Nc == 'Nc6':
    allL = ['12']
    ens_list = 'rt7.5_g0.30'
  elif Nc == 'Nc8':
    allL = ['8', '12', '16']
  else:
    print("ERROR: Unrecognized %s" % Nc)
    sys.exit(1)

  for i in allL:
    L = 'L' + i
    Nt = 'Nt' + i
    vol_dir = Nc + '_' + i + 'nt' + i + '/'
    this_vol = Nc + '/' + Nt + '/' + L
    f.create_group(this_vol)

    # Fourth- and fifth-level groups for each (rt, g) ensemble
    # Set up ensemble lists depending on Nc and L=Nt
    os.chdir(path + vol_dir)
    if Nc == 'Nc4' or Nc == 'Nc6':
      ens_list = ['rt7.5_g0.30']
    else:
      ens_list = glob.glob('rt[5-8]*g0.[123]0*')

    # Add high-temperature Nc8_8nt8 ensembles
    if i == '8':
      for temp_dir in glob.glob('rt[12]*g0.[23]*'):
        ens_list.append(temp_dir)

    for ens in ens_list:
      os.chdir(path + vol_dir + ens)
      toCheck = 'results/poly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished ensembles
        print("Skipping %s" % vol_dir + ens)
        continue

      # Ignore '_c0.0' tag for high-temperature Nc8_8nt8 ensembles
      temp = ens.split('_')
      rt = temp[0]
      zeta = temp[1].replace('g', 'zeta')
      this_ens = this_vol + '/' + rt + '/' + zeta
      this_grp = f.create_group(this_ens)

      # Record t, lambda_lat and mu based on rt and zeta
      Nt_float = float(Nt.split('Nt')[1])
      rt_float = float(rt.split('rt')[1])
      zeta_float = float(zeta.split('zeta')[1])
      this_grp.attrs['t'] = 4.0 / (rt_float * np.sqrt(3.0))
      this_grp.attrs['lambda_lat'] = rt_float / Nt_float
      this_grp.attrs['mu'] = zeta_float * rt_float / Nt_float
      if rt_float > 3.0:
        this_grp.attrs['kappa'] = zeta_float * rt_float / Nt_float
      else:
        this_grp.attrs['kappa'] = 0.0

      # Record thermalization cut and block size
      check = -1
      therm = path + vol_dir + 'therm.sh'
      for line in open(therm):
        if ens in line:
          temp = line.split()
          therm = int(temp[2])
          this_grp.attrs['thermalization cut'] = therm
          this_grp.attrs['block size'] = int(temp[3])
          break
      if therm < 0:
        print("ERROR: Thermalization cut not found for %s: " % this_ens)
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
        this_grp.attrs['autocorrelation time'] = float(temp[0])
      for line in open('results/eig.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time (eig)'] = float(temp[0])
      for line in open('results/bilin.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time (bilin)'] = float(temp[0])

      # Now set up datasets and result-attributes for this ensemble
      # First do plaq, averaging ss & st and also setting Ntraj
      plaq_arr = []
      for line in open('data/plaq.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        dat = 0.5 * (float(temp[1]) + float(temp[2]))
        plaq_arr.append(dat)
      Ntraj = len(plaq_arr)
      this_grp.attrs['Ntraj'] = Ntraj

      plaq = np.array(plaq_arr)
      dset = this_grp.create_dataset('plaq', data=plaq)

      # Record ensemble average as an attribute of the dataset
      for line in open('results/plaq.dat'):
        temp = line.split()
        dset.attrs['ave'] = float(temp[0])
        dset.attrs['err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
          print("%s vs %d in plaq.dat" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Simple observables measured every trajectory=MDTU
      # Bosonic action, link trace, exp(-Delta S)
      # Only want one datum per line from each of these
      # Skipping cg_iters, wall_time
      for obs in ['SB', 'Flink', 'exp_dS']:
        obsfile = 'data/' + obs + '.csv'
        if obs == 'Flink':
          name = 'link_trace'
        else:
          name = obs
        dset = this_grp.create_dataset(name, (Ntraj,), dtype='f')

        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            continue
          dset[traj] = float(temp[1])
          traj += 1

        if traj != Ntraj:
          print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
          print("%d vs %d in %s" % (traj, Ntraj, obsfile))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        for line in open(resfile):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['ave'] = float(temp[0])
          dset.attrs['err'] = float(temp[1])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
            sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Maldacena and Polyakov loops measured every trajectory=MDTU
      # Want three data per line from each of these (real, imag, magnitude)
      # Only magnitude is averaged
      for obs in ['poly_mod', 'poly_mod_polar']:
        obsfile = 'data/' + obs + '.csv'
        if 'polar' in obs:
          name = 'Polyakov_loop'
        else:
          name = 'Maldacena_loop'

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

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        for line in open(resfile):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['magnitude ave'] = float(temp[0])
          dset.attrs['magnitude err'] = float(temp[1])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
            sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Complexified and unitarized Wilson lines measured every trajectory
      # Only monitoring and averaging magnitudes in x- and y-directions
      for obs in ['lines_mod', 'lines_mod_polar']:
        obsfile = 'data/' + obs + '.csv'
        if 'polar' in obs:
          name = 'Unitarized_Wlines'
        else:
          name = 'Complexified_Wlines'

        dset = this_grp.create_dataset(name, (Ntraj,2,), dtype='f')
        dset.attrs['columns'] = ['x-dir magnitude', 'y-dir magnitude']
        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            continue
          dset[traj] = [float(temp[3]), float(temp[4])]
          traj += 1

        if traj != Ntraj:
          print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
          print("%d vs %d in %s" % (traj, Ntraj, obsfile))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        for line in open(resfile):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['x-dir magnitude ave'] = float(temp[2])
          dset.attrs['x-dir magnitude err'] = float(temp[3])
          dset.attrs['y-dir magnitude ave'] = float(temp[4])
          dset.attrs['y-dir magnitude err'] = float(temp[5])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
            sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Violations of fermion bilinear Q-susy Ward identity
      # Save (and count) trajectories at which measurements were run
      bilin_arr = []
      meas_arr = []
      for line in open('data/bilin.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        meas_arr.append(int(temp[0]))
        bilin_arr.append(float(temp[1]))
      Nmeas = len(meas_arr)
      meas = np.array(meas_arr)
      bilin = np.array(bilin_arr)
      this_grp.create_dataset('meas', data=meas)
      dset = this_grp.create_dataset('bilin', data=bilin)

      # Results as attributes, checking Nblocks
      for line in open('results/bilin.dat'):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['ave'] = float(temp[2])
        dset.attrs['err'] = float(temp[3])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
          print("%s vs %d in bilin" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Fermion operator eigenvalues and RHMC spectral range
      # Format: min_eig, max_eig, log(cond_num), min_range, max_range
      # No averaging so attributes are only column labels
      eig = np.empty((Nmeas, 5), dtype = np.float)

      # Minimum eigenvalue from eig.csv
      meas = 0
      for line in open('data/eig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        eig[meas][0] = float(temp[1])
        meas += 1
      if meas != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in min_eig" % (meas, Nmeas))
        sys.exit(1)

      # Maximum eigenvalue from bigeig.csv
      meas = 0
      for line in open('data/bigeig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        eig[meas][1] = float(temp[1])
        meas += 1
      if meas != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in max_eig" % (meas, Nmeas))
        sys.exit(1)

      # Logarithm of condition number
      meas = 0
      for line in open('data/cond_num.csv'):
        if line.startswith('t'):
          continue
        temp = line.split(',')
        eig[meas][2] = float(temp[1])
        meas += 1
      if meas != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in cond_num" % (meas, Nmeas))
        sys.exit(1)

      # Spectral range from order of rational approximation
      meas = 0
      for line in open('data/Norder.csv'):
        if line.startswith('t'):
          continue
        temp = line.split(',')
        srange = spect_range(int(temp[1]))
        eig[meas][3] = float(srange[0])
        eig[meas][4] = float(srange[1])
        meas += 1
      if meas != Nmeas:
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
      # Finally, Wilson line eigenvalue phases
      # Only measured on saved configurations following the thermalization cut
      # Nc phases for all L * Nt sites on z=0 plane
      # !!! Assume all saved configurations separated by 10 MDTU
      Nc_int = int(Nc.split('Nc')[1])
      L_int = int(L.split('L')[1])
      sep = 10
      configs = int((Ntraj - therm) / sep)
      sites = int(L_int * Nt_float)
      WLeig = np.empty((configs, sites, Nc_int), dtype = np.float)

      MDTU = -1
      config = -1
      WLeig_tot = 0
      for line in open('data/WLeig.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        if not int(temp[0]) == MDTU:    # New configuration
          MDTU = int(temp[0])
          config += 1
          site = 0
        for j in range(Nc_int):
          WLeig[config][site][j] = float(temp[j + 1])
          WLeig_tot += 1
        site += 1

      # Too many measurements will cause errors above
      # Also check to make sure no measurements are missing
      expect = configs * sites * Nc_int
      if WLeig_tot != expect:
        print("ERROR: WLeig mismatch in %s: " % this_str, end='')
        print("%d vs %d expected" % (WLeig_tot, expect))
        sys.exit(1)

      dset = this_grp.create_dataset('WLeig_phases', data=WLeig)
      dset.attrs['indices'] = ['config', 'site', 'Nc']
# ------------------------------------------------------------------
