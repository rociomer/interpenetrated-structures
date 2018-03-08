#!/bin/bash

# script for creating a directory for each CSSR structure in the working
#   directory, and submitting the job from within that directory (so as to
#   run on many nodes)

for i in *.cssr;
do
  mkdir ${i%.cssr}

  mv $i ${i%.cssr}/.
  cp *.py ${i%.cssr}/.
  cp interp-generator-submit.sh ${i%.cssr}/.

  cd ${i%.cssr}

  qsub -N ${i%.cssr} interp-generator-submit.sh

  sleep 2m

  cd ../
done
