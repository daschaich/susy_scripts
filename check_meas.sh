#!/bin/bash
counter=0
for dir in /lqcdproj/latticesusy/N4/hybrid_2d/Nc*/rt* ; do
  cd $dir

  # First check that configurations have been generated
  if [ ! -f Out/out.0-10 ] ; then
    continue
  fi

  # Now check whether any measurements are present
  if [ ! -f Out/eig.10 ] ; then
    echo $dir
    continue
  fi

  # Now check that correct number of measurements are present
  o=`ls Out/out.* | wc -l`
  e=`ls Out/eig.* | wc -l`
  c=`ls Out/corr.* | wc -l`
  if [ $o != $e ] || [ $o != $c ] ; then
    echo $dir
  else
    echo "$dir is done"
    let counter++
  fi
done
echo "$counter finished" 
# Alternative: `grep done -c`
