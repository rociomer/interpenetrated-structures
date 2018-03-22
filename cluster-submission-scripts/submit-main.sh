#!/bin/bash

# script for creating a directory for each CSSR structure in the working
#   directory, and submitting the job from within that directory (so as to
#   run on many nodes)

for cssr_file in *.cssr; do
  mkdir ${cssr_file%.cssr}

  mv ${cssr_file} ${cssr_file%.cssr}/.
  cp *.py ${cssr_file%.cssr}/.
  cp interp-generator-submit.sh ${cssr_file%.cssr}/.

  cd ${cssr_file%.cssr}

  qsub -N ${cssr_file%.cssr} interp-generator-submit.sh

  sleep 2m

  cd ../
done
