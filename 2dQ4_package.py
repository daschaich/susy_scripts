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
# Group paths will specify Nc, Nx, Nt,
#   rt = Nt * sqrt(lambda_lat) = 1 / t,
#   and zeta = 'g' = mu / sqrt(lambda_lat) = mu * Nt / rt
# Attributes for each ensemble:
#   t = 1 / rt, lambda_lat, mu, number of trajectories,
#   aspect ratio Nx / Nt, lambda_lat...,
#   mu, number of trajectories, 
#   thermalization cut...,
#   block size...,
#   number of blocks, acceptance rate,
#   autocorrelation times for |ML|, bilinear Ward identity,
#   Tr[X^2] and lowest eigenvalue...
#   Datasets for each ensemble, most with ave and err as attributes:
#   plaquette...,
#   bosonic action, Tr[U.Ubar]/N...,
#   exp(-Delta S), energy density, scalar squares Tr[X^2],
#   real, imag and mag of Maldacena (ML) and Polyakov (PL) loops,
#   real, imag and mag of (complexified and unitarized) spatial Wilson line
#   [with magnitude susceptibility as attribute],
#   fermion bilinear Q-susy Ward identify violation..., [to add]
#   extremal eigenvalues..., [to add]
#   condition number..., [to add]
#   spectral range..., [to add]
#   phases of Wilson line eigenvalues..., [to add]
#   lists of analyzed configurations [to add]
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
      elif Nt == '12':
        allNx = ['12', '24']
      elif Nt == '16':
        allNx = ['12', '24']
      elif Nt == '24':
        allNx = ['24']
      elif Nt == '32':
        allNx = ['32']
      else:
        print("ERROR: Unrecognized Nt=%s for Nc=" % (Nt, Nc))
        sys.exit(1)

    Nt = 'Nt' + Nt
    Nx = 'Nx' + Nx
    vol_dir = Nc + '_' + Nx + 'nt' + Nt + '/'
    this_vol = Nc + '/' + Nt + '/' + Nx
    f.create_group(this_vol)

    # Fourth- and fifth-level groups for each (rt, g) ensemble
    # Set up ensemble lists depending on Nc and L=Nt
    os.chdir(path + vol_dir)
    ens_list = glob.glob('rt*')
    for ens in glob.glob('rt*'):
      temp = ens.split('_')
      rt = temp[0]
      zeta = temp[1]
      this_rtzeta = rt + '/' + zeta
      os.chdir(path + Nc_NxNt + '/' + ens)
      toCheck = 'results/poly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished ensembles
        print("Skipping %s" % Nc_NxNt + '/' + ens)
        continue

    this_ens = this_Nc_NxNt + '/' + this_rtzeta
    this_grp = f.create_group(this_ens)

    # Record lambda_lat and mu_lat based on Nt, t and mu
    rt_float = float(re.split('rt|_g',ens)[1])
    zeta_float = float(ens.split('g')[1])
    Nc_float = float(re.split('Nc|_',Nc_NxNt)[1])
    Nx_float = float(re.split('_|nt',Nc_NxNt)[1])
    Nt_float = float(re.split('_|nt',Nc_NxNt)[2])
    alpha = Nx_float / Nt_float    
    mu = zeta_float * rt_float / Nt_float
    lambda_float = (rt_float / Nt_float)**2  
    temp = 1.0 / rt_float   
    this_grp.attrs['mu'] = mu
    this_grp.attrs['aspect_ratio'] = alpha
    this_grp.attrs['lambda'] = lambda_float
    this_grp.attrs['temperature'] = temp

    # Record thermalization cut and block size
    check = -1
    therm = path + Nc_NxNt + '/therm.sh'
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
      if not line.startswith('W'):
         if not line.startswith('U'):
            temp = line.split()
            this_grp.attrs['autocorrelation time (|ML|)'] = float(temp[0])
    for line in open('results/scalar_sq.autocorr'):
      if not line.startswith('W'):
         if not line.startswith('U'):
            temp = line.split()
            this_grp.attrs['autocorrelation time (Tr[X^2])'] = float(temp[0])

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
    # energy density, plaquette, exp(-Delta S) and Tr[U.Ubar]/N
    # Only want one datum per line from each of these
    for obs in ['energy', 'plaq', 'exp_dS', 'Flink']:
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
    # Unitarized Wilson line also measured every trajectory=MDTU
    # Want three data per line from each of these (real, imag, magnitude)
    # Only magnitude is averaged
    for obs in ['lines_mod_polar']:
      obsfile = 'data/' + obs + '.csv'

      dset = this_grp.create_dataset(obs, (Ntraj,3,), dtype='f')
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
    # Maldacena loop, Polyakov loop, and Wilson line also measured every trajectory=MDTU
    # Want three data per line from each of these (real, imag, magnitude)
    # Only magnitude is averaged
    for obs in ['poly_mod', 'poly_mod_polar', 'lines_mod']:
      obsfile = 'data/' + obs + '.csv'

      dset = this_grp.create_dataset(obs, (Ntraj,3,), dtype='f')
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
    # Want four data per line (extremal minimum and maximum eigenvalues)
    for obs in ['scalar_eig']:
      obsfile = 'data/' + obs + '.csv'

      dset = this_grp.create_dataset(obs, (Ntraj,4,), dtype='f')
      dset.attrs['columns'] = ['min_min', 'min_max', 'max_min', 'max_max']
      traj = 0
      for line in open(obsfile):
        temp = line.split(',')
        if line.startswith('M') or line.startswith('t'):
          continue
        dset[traj] = [float(temp[1]), float(temp[2]), float(temp[3]), float(temp[4])]
        traj += 1

      if traj != Ntraj:
        print("ERROR: Ntraj mismatch in %s, " % this_ens, end='')
        print("%d vs %d in %s" % (traj, Ntraj, obsfile))
        sys.exit(1)
    # ------------------------------------------------------------
  
    # ------------------------------------------------------------
    # Now observables which are observed every 10 MDTU
    # First do scalar squares to set new Ntraj
    scalarsq_arr = []
    for line in open('data/scalar_sq.csv'):
      if line.startswith('M'):
        continue
      temp = line.split(',')
      scalarsq_arr.append(float(temp[1]))
    Ntraj = len(scalarsq_arr)
    this_grp.attrs['Ntraj'] = Ntraj

    scalarsq = np.array(scalarsq_arr)
    dset = this_grp.create_dataset('scalarsq', data=scalarsq)

    # Record ensemble average as an attribute of the dataset
    for line in open('results/scalar_sq.dat'):
      temp = line.split()
      dset.attrs['ave'] = float(temp[0])
      dset.attrs['err'] = float(temp[1])
      if not int(temp[-1]) == Nblocks:
        print("ERROR: Nblocks mismatch in %s: " % this_ens, end='')
        print("%s vs %d in SB.dat" % (temp[-1], Nblocks))
        sys.exit(1)
    # ------------------------------------------------------------
