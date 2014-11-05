#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# Parse N=(2, 2) SYM output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting
# Normalize the Polyakov loop data by Nc and the bosonic action by 3Nc^2 / 2

die "Usage: $0 <path>\n"
  if (@ARGV != 1);

my $dir = shift;
my $path = "/nfs/beowulf02/schaich/SYM/$dir";
my $i;
my $junk;
my $temp;
open ERRFILE, ">> ERRORS" or die "Error opening ERRORS ($!)\n";
open MISSINGFILES, "> MISSING" or die "Error opening MISSING ($!)\n";

# Physical observables
open PLAQ, "> $path/data/plaq.csv" or die "Error opening $path/data/plaq.csv ($!)\n";
print PLAQ "MDTU,plaq\n";
open SB, "> $path/data/SB.csv" or die "Error opening $path/data/SB.csv ($!)\n";
print SB "MDTU,S_B\n";
open POLY, "> $path/data/poly.csv" or die "Error opening $path/data/poly.csv ($!)\n";
print POLY "ReTr(L),ImTr(L)\n";
open POLY_MOD, "> $path/data/poly_mod.csv" or die "Error opening $path/data/poly_mod.csv ($!)\n";
print POLY_MOD "MDTU,|Tr(L)|,ReTr(L),ImTr(L)\n";
open FLINK, "> $path/data/Flink.csv" or die "Error opening $path/data/Flink.csv ($!)\n";
print FLINK "MDTU,link\n";
open DET, "> $path/data/det.csv" or die "Error opening $path/data/det.csv ($!)\n";
print DET "MDTU,|det-1|,|Re(det)-1|\n";
open EIG, "> $path/data/eig.csv" or die "Error opening $path/data/eig.csv ($!)\n";
print EIG "MDTU,0,2,4,6,8,10\n";
open BILIN, "> $path/data/bilin.csv" or die "Error opening $path/data/bilin.csv ($!)\n";
print BILIN "MDTU,susyTrans\n";

# Evolution observables
open ACCP, "> $path/data/accP.csv" or die "Error opening $path/data/accP.csv ($!)\n";
print ACCP "t,accP\n";
open EXP_DS, "> $path/data/exp_dS.csv" or die "Error opening $path/data/exp_dS.csv ($!)\n";
print EXP_DS "t,e^(-dS)\n";
open DELTAS, "> $path/data/deltaS.csv" or die "Error opening $path/data/deltaS.csv ($!)\n";
print DELTAS "t,deltaS\n";
open ABS_DS, "> $path/data/abs_dS.csv" or die "Error opening $path/data/abs_dS.csv ($!)\n";
print ABS_DS "t,|deltaS|\n";
open FORCE, "> $path/data/force.csv" or die "Error opening $path/data/force.csv ($!)\n";
print FORCE "t,G,F\n";
open CG_ITERS, "> $path/data/cg_iters.csv" or die "Error opening $path/data/cg_iters.csv ($!)\n";
print CG_ITERS "t,cg_iters\n";
open WALLTIME, "> $path/data/walltime.csv" or die "Error opening $path/data/walltime.csv ($!)\n";
print WALLTIME "t,walltime\n";
open WALLTU, "> $path/data/wallTU.csv" or die "Error opening $path/data/wallTU.csv ($!)\n";
print WALLTU "t,cost\n";

# Run parameters
open NSTEP, "> $path/data/Nstep.csv" or die "Error opening $path/data/Nstep.csv ($!)\n";
print NSTEP "t,N_f,N_g\n";
open STEPSIZE, "> $path/data/stepsize.csv" or die "Error opening $path/data/stepsize.csv ($!)\n";
print STEPSIZE "t,eps_f,eps_g\n";
open TLENGTH, "> $path/data/tlength.csv" or die "Error opening $path/data/tlength.csv ($!)\n";
print TLENGTH "t,L\n";
open KEY, "> $path/data/key.csv" or die "Error opening $path/data/key.csv ($!)\n";
print KEY "t,file\n";
open TU, "> $path/data/TU.csv" or die "Error opening $path/data/TU.csv ($!)\n";
print TU "t,MDTU\n";
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Declare observables to be printed
my $cpus;
my $stepsize;
my $stepsize_gauge;
my $Nstep;
my $Nstep_gauge;
my $tlength;
my $acc;
my $force_g;
my $force_f;
my $dS;
my $exp_dS;
my $abs_dS;
my $iters;
my $plaq;
my $b_act;
my $ploop_r;
my $ploop_i;
my $p_mod;
my $ave_link;
my $det_r;
my $det_i;
my $det;
my $detSq_r;
my $detSq_i;
my $ave_time;
my $TUtime;
my $bilin;
my $trace;
my $gauge;
my $susy;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# These fellows help track status and perform checks
my $infile = -1;
my $check = -1;           # Check that file is present or completed
my $walltime = -1;        # Check that file completed successfully
my $load;                 # MDTU of loaded configuration
my $cfg;                  # MDTU of saved configuration
my $oldcfg = 0;           # Check if any files are missing
my $stamp = "start";
my $meas_stamp = "";
my $oldstamp = "start";   # Check that correct configuration was used
my $starting;             # Check if this is the start of the job
my $Nc = -1;              # For Polyakov loop normalization

# Running sums for the ensemble as a whole
my $traj = 0;
my $endtraj = 0;          # Reset from traj, so not really running sum
my $MDTU = 0;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through files out.$load-$save from list.txt
open FILES, "< list.txt" or die "Error opening list.txt ($!)\n";
my @files = <FILES>;
close FILES;
FILE: for my $file (@files) {
  chomp ($file);          # Get rid of linebreaks in filenames
  ($load, $cfg) = split /-/, $file;
  # Initialize running sums and set dummy walltime
  # If the walltime isn't overwritten, then the run died
  # or its output file is corrupted
  $walltime = -1;
  $stamp = "start";

  # Open file
  # If not found, move on to next file instead of killing whole program,
  # but print error message so I know there is a problem
  $infile = "$path/Out/out.$file";
  $check = open IN, "< $infile";
  if (!$check) {
    print STDERR "Problem opening $path/Out/out.$file: $!\n";
    print ERRFILE "Problem opening $path/Out/out.$file: $!\n";
    next FILE;
  }
#  print STDOUT "$infile\n"; # Monitor running status
  my @in = <IN>;
  close IN;

  # If not starting from file 1 in this ensemble,
  # or if we seem to have skipped a file,
  # guess approximate starting trajectory
  my $traj_per_file = -1;
  LINE: for my $line (@in) {
    if ($line =~/^PLACEHOLDER/) {
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      $walltime = -2;
    }
    # Extract Nc for bosonic action and Polyakov loop normalizations
    elsif ($line =~ /^N=\(2, 2\) SYM, /) {
      ($junk, $junk, $temp, $junk) = split /,/, $line;
      # Funny business due to leading space
      ($junk, $junk, $junk, $Nc) = split /\s+/, $temp;
    }
    elsif ($line =~ /^trajecs /) {
      ($junk, $traj_per_file) = split /\s+/, $line;
      $endtraj = $traj + $traj_per_file;
      last LINE;   # Don't go through whole file yet
    }
  }
  if ($traj_per_file == -1) {
    print STDERR "$infile never defines number of trajectories\n";
    print ERRFILE "$infile never defines number of trajectories\n";
    next FILE;     # Skip this file, which appears malformed
  }

  if ($traj == 0 && $load > 0 || $load != $oldcfg) {
    print STDERR "$infile misses MDTU before $load\n";
    $traj = $load;
    $MDTU = $load;
    $endtraj = $traj + $traj_per_file;
  }
  # --------------------------------------------------------------



  # --------------------------------------------------------------
  # Cycle through lines in the file
  for my $line (@in) {
    # Check for premature termination (e.g., layout problem)
    if ($line =~ /^termination/) {
      print STDERR "$infile reports premature termination\n";
      print ERRFILE "$infile reports premature termination\n";
      next FILE;     # Skip this file, which appears malformed
    }
  }

  # At this point, we should be able to begin
  $oldcfg = $cfg;
  for my $line (@in) {
    # Extract constant run parameters
    if ($line =~ /^traj_length /) {
      ($junk, $tlength) = split /\s+/, $line;
      print TLENGTH "$endtraj,$tlength\n";
    }

    elsif ($line =~ /^nstep /) {
      ($junk, $Nstep) = split /\s+/, $line;
      $stepsize = $tlength / $Nstep;
    }
    elsif ($line =~ /^nstep_gauge /) {
      ($junk, $Nstep_gauge) = split /\s+/, $line;
      $stepsize_gauge = $stepsize / (2 * $Nstep_gauge);
      $Nstep_gauge *= 2 * $Nstep;
      print NSTEP "$endtraj,$Nstep,$Nstep_gauge\n";
      print STEPSIZE "$endtraj,$stepsize,$stepsize_gauge\n";
    }

    elsif ($line =~ /^Machine = /) {
      ($junk, $junk, $junk, $junk, $junk, $cpus, $junk) = split /\s+/, $line;
    }

    # Check that the file loaded the appropriate configuration
    elsif ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      if ($stamp eq "start") {    # Loading configuration
        $stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
        if ($stamp ne $oldstamp && $oldstamp ne "start") {
          print STDERR "$infile: loaded $stamp doesn't match saved $oldstamp\n";
          print ERRFILE "$infile: loaded $stamp doesn't match saved $oldstamp\n";
        }
      }
      else {                # Saving configuration
        $oldstamp = join ' ', split ' ', $line;  # Prepare for next file or MCRG file
        $stamp = "start";
      }
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Now extract evolution observables and physical observables
    # Acceptance comes before measurements
    elsif ($line =~ / delta S = /) {
      $traj++;            # Remains across files in this ensemble!
      $MDTU += $tlength;  # Remains across files in this ensemble!
      # Round off digits in MDTU
      $junk = $MDTU;
      $MDTU = sprintf("%.3f", $junk);
      print TU "$traj,$MDTU\n";
      print KEY "$MDTU,$file\n";

      ($acc, $junk, $junk, $junk, $temp, $junk) = split /\s+/, $line;
      # Strip nasty comma from dS!
      if ($temp =~ /,/) {
        ($dS, $junk) = split /,/, $temp;
      }
      else {
        $dS = $temp
      }
      print DELTAS "$traj,$dS\n";
      $exp_dS = exp(-$dS);
      print EXP_DS "$traj,$exp_dS\n";

      # For RMS, don't have an easy way to average over many measurements
      # Instead just print out absolute value and consider its running average
      $abs_dS = abs($dS);
      print ABS_DS "$traj,$abs_dS\n";

      # Will be smeared out by running averages
      if ($acc =~ /ACCEPT/) {
        print ACCP "$traj,1\n";
      }
      else {
        print ACCP "$traj,0\n";
      }
    }

    # Forces
    elsif ($line =~ /MONITOR_FORCE_GAUGE /) {
      ($junk, $force_g) = split /\s+/, $line;
    }
    elsif ($line =~ /MONITOR_FORCE_FERMION0 /) {
      ($junk, $force_f) = split /\s+/, $line;
      print FORCE "$traj,$force_g,$force_f\n";
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Average value of the link -- ignore the individual directions for now
    # Printed twice at start of job -- ignore one
    elsif ($line =~ /^START /) {
      $starting = 1;
    }
    elsif ($line =~ /^FLINK /) {
      if ($starting == 1) {
        $starting = 0;
      }
      else {
        ($junk, $junk, $junk, $ave_link) = split /\s+/, $line;
        print FLINK "$MDTU,$ave_link\n";
      }
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Gauge measurements come next:
    # plaquette, Polyakov loop, bosonic action, and CG iterations for some reason
    elsif ($line =~ /^GMES /) {
      ($junk, $ploop_r, $ploop_i, $iters, $plaq, $b_act) = split /\s+/, $line;
      print PLAQ "$MDTU,$plaq\n";
      print CG_ITERS "$traj,$iters\n";

      # Normalize bosonic action and Polyakov loop using Nc extracted above
      $b_act /= 1.5 * $Nc**2;
      print SB "$MDTU,$b_act\n";

      $ploop_r /= $Nc;
      $ploop_i /= $Nc;
      print POLY "$ploop_r,$ploop_i\n";
      $p_mod = sqrt($ploop_r**2 + $ploop_i**2);
      print POLY_MOD "$MDTU,$p_mod,$ploop_r,$ploop_i\n";
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Finally, the determinant is only measured every once in a while
    elsif ($line =~ /^DET /) {
      ($junk, $temp, $det_i, $detSq_r, $detSq_i) = split /\s+/, $line;
      $det_r = abs($temp - 1.0);
      $det = sqrt($det_r**2 + $det_i**2);
      print DET "$MDTU,$det,$det_r\n";
    }

    # Store total walltime to average at the end
    elsif ($line =~ /^Time = /) {
      ($junk, $junk, $walltime, $junk) = split /\s+/, $line;
    }
  } # Done cycling through output file
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Check to see if run seems to have finished properly
  if ($walltime == -1) {
    print STDERR "$infile didn't print final timing\n";
    print ERRFILE "$infile didn't print final timing\n";
  }
  elsif ($walltime == -2) {
    # Placeholder file -- error has been addressed as well as possible,
    # but don't print nonsense wall clock time
  }
  else {    # We are good to go
    # Average walltime over all trajectories
    $ave_time = $walltime / $traj_per_file;
    print WALLTIME "$traj,$ave_time\n";
    $TUtime = $ave_time * $cpus / (60 * abs($tlength));
    print WALLTU "$traj,$TUtime\n";
  }
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
EIGEN:
  # Now deal with the corresponding eigenvalues file, if it is there
  $check = -1;
  $infile = "$path/Out/eig.$cfg";
  $check = open EIG_IN, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
    print MISSINGFILES "$infile\n";
    goto CORR;
  }
  my @eig_in = <EIG_IN>;
  close EIG_IN;

  # We have a file, so let's cycle over its lines
  # These are always paired (checked by check_eig_pairs.py)
  # Focus on first six pairs, 0, 2, 4, 6, 8 and 10
  $check = -1;               # Check whether file completed successfully
  my @eig = ("null", "null", "null", "null", "null", "null");
  for my $line (@eig_in) {
    if ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      $meas_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
      if ($meas_stamp ne $oldstamp) {
        print STDERR "$infile: loaded $meas_stamp doesn't match saved $oldstamp\n";
        print ERRFILE "$infile: loaded $meas_stamp doesn't match saved $oldstamp\n";
      }
    }
    elsif ($line =~ /^EIGENVALUE 0 /) {
      ($junk, $junk, $eig[0], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 2 /) {
      ($junk, $junk, $eig[1], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 4 /) {
      ($junk, $junk, $eig[2], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 6 /) {
      ($junk, $junk, $eig[3], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 8 /) {
      ($junk, $junk, $eig[4], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 10 /) {
      ($junk, $junk, $eig[5], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /WARNING/) {
      print STDERR "$infile saturated eigenvalue iterations\n";
      print ERRFILE "$infile saturated eigenvalue iterations\n";
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done with eigenvalues file
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }
  print EIG "$MDTU,$eig[0],$eig[1],$eig[2],$eig[3],$eig[4],$eig[5]\n";
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
CORR:
  # Finally deal with any other corresponding measurements
  # For now just fermion bilinears and related quantities
  $check = -1;
  $infile = "$path/Out/corr.$cfg";
  $check = open CORR_IN, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
    print MISSINGFILES "$infile\n";
    next FILE;
  }
  my @corr_in = <CORR_IN>;
  close CORR_IN;

  # Need to add factor of C2 to gauge piece
  # Try to extract C2 from $path
  my $C2 = 1.0;
  if ($path =~ /-c/) {
    ($junk, $C2) = split /-c/, $path;
  }

  # We have a file, so let's cycle over its lines
  $check = -1;               # Check whether file completed successfully
  for my $line (@corr_in) {
    if ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      $meas_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
      if ($meas_stamp ne $oldstamp) {
        print STDERR "$infile: loaded $meas_stamp doesn't match saved $oldstamp\n";
        print ERRFILE "$infile: loaded $meas_stamp doesn't match saved $oldstamp\n";
      }
    }
#    elsif ($line =~ /^BILIN /) {
#      ($junk, $bilin, $junk) = split /\s+/, $line;
#    }
    elsif ($line =~ /^SUSY /) {
      ($junk, $trace, $junk, $gauge, $temp) = split /\s+/, $line;
      $susy = ($C2 * $gauge - $trace) / ($C2 * $gauge + $trace);
#      $trace -= $bilin;                     # Isolate U(1) piece
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done with measurements file
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }
#  print BILIN "$MDTU,$bilin,$trace,$gauge\n";
  print BILIN "$MDTU,$susy\n";
} # Done cycling through files
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
close ERRFILE;
close MISSINGFILES;
close KEY;
close STEPSIZE;
close NSTEP;
close TLENGTH;
close TU;
close DELTAS;
close EXP_DS;
close ABS_DS;
close ACCP;
close FORCE;
close PLAQ;
close SB;
close POLY;
close POLY_MOD;
close FLINK;
close DET;
close CG_ITERS;
close WALLTIME;
close WALLTU;
close EIG;
close BILIN;
exit(0);
# ------------------------------------------------------------------
