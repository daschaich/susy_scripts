#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package bosonic BMN data and results/attributes into HDF5 file

# Cycle over ensembles and write to ~/SYM/bosonicBMN/bBMN_data.h5
# A total of 397 ensembles, all with Nt=24
# Group paths will specify Nc, mu = mu_lat / lambda_lat^(1/3)
#   and t = 1 / (Nt * lambda_lat^(1/3))
# Attributes for each ensemble:
#   lambda_lat, mu_lat, number of trajectories, thermalization cut,
#   block size, number of blocks, acceptance rate,
#   autocorrelation times for |PL| and one Tr[X^2]
# Datasets for each ensemble, most with ave and err as attributes:
#   bosonic action, exp(-Delta S), scalar squares Tr[X^2],
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
for Nc_mu in glob.glob('Nc*'):
  temp = Nc_mu.split('_24_')
  Nc = temp[0]
  mu = temp[1]
  this_Nc_mu = Nc + '/' + mu
  f.create_group(this_Nc_mu)

  # Third-level group for each ensemble specified by t
  os.chdir(path + Nc_mu)
  for ens in glob.glob('t*'):
    os.chdir(path + Nc_mu + '/' + ens)
    toCheck = 'results/poly_mod.autocorr'
    if not os.path.isfile(toCheck):     # Skip unfinished ensembles
      print("Skipping %s" % Nf + '/' + vol + '/' + ens)
      continue

    temp = ens.split('_')
    g = temp[0]
    f = temp[1]
    this_str = this_Nc_mu + '/' + g + '/' + f
    this_grp = f.create_group(this_str)

    # Record thermalization cut and block size
    # TODO: Check Out/analysis.sh Out/observables.sh ...
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
      print("ERROR: Thermalization cut not found for %s: " % this_str)
      sys.exit(1)

    # Record acceptance, auto-correlation and Nblocks for ensemble
    for line in open('results/accP.dat'):
      temp = line.split()
      Nblocks = int(temp[-1])
      this_grp.attrs['Nblocks'] = Nblocks
      this_grp.attrs['acceptance'] = float(temp[0])
    for line in open('results/poly_mod.autocorr'):
      temp = line.split()
      this_grp.attrs['autocorrelation time'] = float(temp[0])
    for line in open('results/scalarsquares.autocorr'):
      if not line.startswith('W'):      # Ignore warning message
        temp = line.split()
        this_grp.attrs['autocorrelation time (Tr[x^2])'] = float(temp[0])

# TODO...
      # Now set up data and attributes/results for this ensemble
      # First do simple observables measured every trajectory=MDTU
      # All these have a single datum per line
      # Skipping wall_time
      for obs in ['poly_mod', 'poly_arg', 'poly_r', 'exp_dS']:
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

        # Same for susceptibility, skewness, kurtosis if present
        suscfile = 'results/' + obs + '.suscept'
        if os.path.isfile(suscfile):        # No suscept for some obs
          for line in open(suscfile):
            temp = line.split()
            dset.attrs[temp[0]] = float(temp[1])
            dset.attrs[temp[0] + '_err'] = float(temp[2])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
              print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
              sys.exit(1)
# ------------------------------------------------------------------
