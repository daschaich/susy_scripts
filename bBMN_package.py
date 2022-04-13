#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package bosonic BMN data and results/attributes into HDF5 file

# Cycle over ensembles and write to ~/SYM/bosonicBMN/bBMN_data.h5
# A total of 384 ensembles, all with Nt=24
# Group paths will specify Nc, mu = mu_lat / lambda_lat^(1/3)
#   and t = 1 / (Nt * lambda_lat^(1/3))
# Attributes for each ensemble:
#   lambda_lat, mu_lat, number of trajectories, thermalization cut,
#   block size, number of blocks, acceptance rate,
#   autocorrelation times for |PL| and one Tr[X^2]
# Datasets for each ensemble, most with ave and err as attributes:
#   bosonic action, exp(-Delta S), Myers term, scalar squares Tr[X^2],
#   real, imag and magnitude of Polyakov loop (PL),
#   phases of PL eigenvalues, lists of analyzed configurations
# PL magnitude susceptibility also included as dataset attribute
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Annoyingly, need to use absolute Barkla path
f = h5py.File("/users/schaich/SYM/bosonicBMN/bBMN_data.h5", 'w')
path = "/users/schaich/anosh/susy_bound/1d/b_g_bmn/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)

# Top- and second-level groups for each Nc and mu
Nt = 24.0
for Nc_mu in glob.glob('Nc*'):
  temp = Nc_mu.split('_24_')
  Nc = temp[0]
  mu = temp[1]
  this_Nc_mu = Nc + '/' + mu
  f.create_group(this_Nc_mu)

  # Third-level group for each ensemble specified by t
  os.chdir(path + Nc_mu)
  for ens in glob.glob('t[0-9]*'):
    os.chdir(path + Nc_mu + '/' + ens)
    toCheck = 'results/poly_mod.autocorr'
    if not os.path.isfile(toCheck):     # Skip unfinished ensembles
      print("Skipping %s" % Nc_mu + '/' + ens)
      continue

    this_ens = this_Nc_mu + '/' + ens
    this_grp = f.create_group(this_ens)

    # Record lambda_lat and mu_lat based on Nt, t and mu
    t_float = float(ens.split('t')[1])
    mu_float = float(mu.split('mu')[1])
    cube_root = 1.0 / (Nt * t_float)            # lambda_lat^{1/3}
    this_grp.attrs['lambda_lat'] = np.power(cube_root, 3)
    this_grp.attrs['mu_lat'] = mu_float * cube_root

    # Record thermalization cut and block size
    check = -1
    therm = path + Nc_mu + '/therm.sh'
    for line in open(therm):
      if ens in line:
        temp = line.split()
        this_grp.attrs['thermalization cut'] = int(temp[2])
        this_grp.attrs['block size'] = int(temp[3])
        check = int(temp[2])
        break
    if check < 0:
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
      if not line.startswith('W'):      # Ignore warning message
        temp = line.split()
        this_grp.attrs['autocorrelation time'] = float(temp[0])
    for line in open('results/scalarsquares.autocorr'):
      if not line.startswith('W'):      # Ignore warning message
        temp = line.split()
        this_grp.attrs['autocorrelation time (Tr[x^2])'] = float(temp[0])

    # Now set up datasets and result-attributes for this ensemble
    # First do bosonic action to set Ntraj
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
        print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
        print("%s vs %d in SB.dat" % (temp[-1], Nblocks))
        sys.exit(1)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Simple observables measured every trajectory=MDTU
    # SO(6)/SO(3) ratio, Myers term, extent of space, exp(-Delta S)
    # Only want one datum per line from each of these
    # Myers and extent of space normalized while parsing
    for obs in ['ratio', 'Myers', 'extent', 'exp_dS']:
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
        print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
        print("%d vs %d in %s" % (traj, Ntraj, obsfile))
        sys.exit(1)

      # Results as attributes, checking Nblocks
      resfile = 'results/' + obs + '.dat'
      if os.path.isfile(resfile):
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
    # Internal energy also measured every trajectory=MDTU
    # Want two data per line: E/N^2 and E_prime (normalized while parsing)
    # Only internal energy is averaged
    for obs in ['poly_mod']:
      obsfile = 'data/' + obs + '.csv'
      name = 'Polyakov_loop'

      dset = this_grp.create_dataset(name, (Ntraj,2,), dtype='f')
      dset.attrs['columns'] = ['E/N^2', 'E_prime']
      traj = 0
      for line in open(obsfile):
        temp = line.split(',')
        if line.startswith('M') or line.startswith('t'):
          continue
        dset[traj] = [float(temp[1]), float(temp[2])]
        traj += 1

      if traj != Ntraj:
        print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
        print("%d vs %d in %s" % (traj, Ntraj, obsfile))
        sys.exit(1)

      # Results (including specific heat) as attributes, checking Nblocks
      resfile = 'results/' + obs + '.dat'
      for line in open(resfile):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['energy ave'] = float(temp[0])
        dset.attrs['energy err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
          print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
          sys.exit(1)
      specfile = 'results/' + obs + '.specheat'
      for line in open(specfile):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['specific heat'] = float(temp[0])
        dset.attrs['specific heat_err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
          print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
          sys.exit(1)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Polyakov loop also measured every trajectory=MDTU
    # Want three data per line from each of these (real, imag, magnitude)
    # Only magnitude is averaged
    for obs in ['poly_mod']:
      obsfile = 'data/' + obs + '.csv'
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
      suscfile = 'results/' + obs + '.suscept'
      for line in open(suscfile):
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

    # Too many measurements will cause errors above
    # Also check to make sure no measurements are missing
    expect = Ntraj * Nscalar
    if TrXSq_tot != expect:
      print("ERROR: TrXSq mismatch in %s: " % this_ens, end='')
      print("%d vs %d expected" % (TrXSq_tot, expect))
      sys.exit(1)

    dset = this_grp.create_dataset('Scalar_squares', data=TrXSq)

    # Results as attributes, checking Nblocks
    resfile = 'results/scalarsquares.dat'
    for line in open(resfile):
      TrXSq_ave = np.empty((Nscalar), dtype = np.float)
      TrXSq_err = np.empty_like(TrXSq_ave)
      temp = line.split()
      for j in range(Nscalar):
        TrXSq_ave[j] = float(temp[2 * j])
        TrXSq_err[j] = float(temp[2 * j + 1])

      dset.attrs['magnitude ave'] = TrXSq_ave
      dset.attrs['magnitude err'] = TrXSq_err
      if not int(temp[-1]) == Nblocks:
        print("ERROR: Nblocks mismatch in %s, " % this_ens, end='')
        print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
        sys.exit(1)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Finally, Polyakov loop eigenvalue phases
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

    # Too many measurements will cause errors above
    # Also check to make sure no measurements are missing
    expect = Ntraj * Nc_int
    if PLeig_tot != expect:
      print("ERROR: PLeig mismatch in %s: " % this_ens, end='')
      print("%d vs %d expected" % (PLeig_tot, expect))
      sys.exit(1)

    dset = this_grp.create_dataset('PLeig_phases', data=PLeig)
# ------------------------------------------------------------------
