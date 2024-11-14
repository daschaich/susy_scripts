#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package 2d N=(2,2) SYM data and results/attributes into HDF5 file

# Cycle over ensembles and write to /mnt/lustre/users/anosh/susy_bound/susy_bound/2dQ4_data.h5
# 207 Nc=12 ensembles plus 12 Nc=16 12nt12 and 12 Nc=20 12nt12 ensembles
# Group paths will specify Nc, Nx, Nt, rt = Nt * sqrt(lambda_lat),
#   and zeta = 'g' = mu / sqrt(lambda_lat) = mu * Nt / rt
# Attributes for each ensemble:
#   aspect ratio Nx / Nt, t = 1 / rt, lambda_lat, mu, number of trajectories,
#   thermalization cut, block size, number of blocks, acceptance rate,
#   and autocorrelation times for |ML|, lowest fermion operator eigenvalue,
#     bilinear Ward identity and Tr[X^2]
# Datasets for each ensemble, most with ave and err as attributes:
#   plaquette, bosonic action, energy density, Tr[U.Ubar]/N, exp(-Delta S),
#   mag of (complexified and unitarized) spatial Wilson line
#     (with unitarized susceptibility as attribute),
#   real, imag and mag of Maldacena (ML) and Polyakov (PL) loops,
#   extremal scalar eigenvalues,
#   fermion bilinear Q-susy Ward identify violation,
#   list of analyzed configurations,
#   scalar square Tr[X^2], extremal fermion operator eigenvalues,
#   condition number and spectral range
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
path = "/mnt/lustre/users/anosh/susy_bound/susy_bound/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)
f = h5py.File("/mnt/lustre/users/anosh/susy_bound/susy_bound/2dQ4_data.h5", 'w')

# Top-level groups for each Nc
for Nc in ['Nc12', 'Nc16', 'Nc20']:
  f.create_group(Nc)

  # Second- and third-level groups for each Nt and Nx
  # Currently only one volume for Nc=16 and 20, but let's be a bit flexible
  if Nc == 'Nc16' or Nc == 'Nc20':
    allNt = ['12']
  elif Nc == 'Nc12':
    allNt = ['6', '12', '16', '24', '32']
  else:
    print("ERROR: Unrecognized Nc=%s" % Nc)
    sys.exit(1)

  for Nt in allNt:
    if Nc == 'Nc16' or Nc == 'Nc20':
      allNx = ['12']
    elif Nc == 'Nc12':
      if Nt == '6':
        allNx = ['24']
      elif Nt == '12' or Nt == '24':
        allNx = ['12', '24']
      elif Nt == '16':
        allNx = ['16', '24']
      elif Nt == '32':
        allNx = ['32']
      else:
        print("ERROR: Unrecognized Nt=%s for Nc=" % (Nt, Nc))
        sys.exit(1)

    Nt_tag = 'Nt' + Nt
    for Nx in allNx:
      Nx_tag = 'Nx' + Nx
      vol_dir = Nc + '_' + Nx + 'nt' + Nt + '/'
      this_vol = Nc + '/' + Nt_tag + '/' + Nx_tag
      f.create_group(this_vol)

      # Fourth- and fifth-level groups for each (rt, g) ensemble
      # Set up ensemble lists for each Nc and volume
      os.chdir(path + vol_dir)
      for ens in glob.glob('rt*'):
        os.chdir(path + vol_dir + ens)
        toCheck = 'results/poly_mod.autocorr'
        if not os.path.isfile(toCheck):     # Skip unfinished ensembles
          print("Skipping %s" % vol_dir + ens)
          continue

        temp = ens.split('_')
        rt = temp[0]
        zeta = temp[1].replace('g', 'zeta')
        this_ens = this_vol + '/' + rt + '/' + zeta
        this_grp = f.create_group(this_ens)

        # Record alpha, t, lambda_lat and mu based on Nt, Nx, rt and zeta
        Nt_float = float(Nt_tag.split('Nt')[1])
        Nx_float = float(Nx_tag.split('Nx')[1])
        rt_float = float(rt.split('rt')[1])
        zeta_float = float(zeta.split('zeta')[1])
        this_grp.attrs['alpha'] = Nx_float / Nt_float
        this_grp.attrs['t'] = 1.0 / rt_float
        this_grp.attrs['lambda_lat'] = (rt_float / Nt_float)**2
        this_grp.attrs['mu'] = zeta_float * rt_float / Nt_float

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
          this_grp.attrs['autocorrelation time (|ML|)'] = float(temp[0])
        for line in open('results/eig.autocorr'):
          temp = line.split()
          this_grp.attrs['autocorrelation time (eig)'] = float(temp[0])
        for line in open('results/bilin.autocorr'):
          temp = line.split()
          this_grp.attrs['autocorrelation time (bilin)'] = float(temp[0])
        for line in open('results/scalar_sq.autocorr'):
          temp = line.split()
          this_grp.attrs['autocorrelation time (Tr[X^2])'] = float(temp[0])

        # Now set up datasets and result-attributes for this ensemble
        # First do plaquette and set Ntraj
        plaq_arr = []
        for line in open('data/plaq.csv'):
          if line.startswith('M'):
            continue
          temp = line.split(',')
          plaq_arr.append(float(temp[1]))
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
        # Bosonic action, energy density, Tr[U.Ubar]/N, exp(-Delta S),
        # and complexified and unitarized Wilson line magnitudes
        # Only want one datum per line from each of these
        # Skipping cg_iters, wall_time
        for obs in ['SB', 'energy', 'Flink', 'exp_dS', \
                    'lines_mod', 'lines_mod_polar']:
          obsfile = 'data/' + obs + '.csv'
          if obs == 'Flink':
            name = 'link_trace'
          elif obs == 'lines_mod':
            name = 'Complexified_Wlines'
          elif obs == 'lines_mod_polar':
            name = 'Unitarized_Wlines'
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

          # Include unitarized Wilson line susceptibility
          if obs == 'lines_mod_polar':
            suscfile = 'results/' + obs + '.suscept'
            for line in open(suscfile):
              if line.startswith('#'):
                continue
              temp = line.split()
              dset.attrs['suscept'] = float(temp[0])
              dset.attrs['suscept_err'] = float(temp[1])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
              print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
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
        # Extremal scalar eigenvalues also measured every trajectory=MDTU
        # Want two data per line (average minimum and maximum eigenvalues)
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

        # Results as attributes, checking Nblocks
        for line in open('results/scalar_eig_ave.dat'):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['min ave'] = float(temp[0])
          dset.attrs['min err'] = float(temp[1])
          dset.attrs['max ave'] = float(temp[2])
          dset.attrs['max err'] = float(temp[3])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in scalar_eig_ave" % (temp[-1], Nblocks))
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
          dset.attrs['ave'] = float(temp[0])
          dset.attrs['err'] = float(temp[1])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in bilin" % (temp[-1], Nblocks))
            sys.exit(1)
        # ------------------------------------------------------------

        # ------------------------------------------------------------
        # Next do average Tr[X^2]
        dset = this_grp.create_dataset('scalar_sq', (Nmeas,), dtype='f')
        meas = 0
        for line in open('data/scalar_sq.csv'):
          temp = line.split(',')
          if line.startswith('M'):
            continue
          dset[meas] = float(temp[1])
          meas += 1

        if meas != Nmeas:
          print("ERROR: Nmeas mismatch in %s, " % this_ens, end='')
          print("%d vs %d in scalar_sq" % (meas, Nmeas))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        for line in open('results/scalar_sq.dat'):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['ave'] = float(temp[0])
          dset.attrs['err'] = float(temp[1])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
            print("%s vs %d in scalar_sq" % (temp[-1], Nblocks))
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
