#!/bin/bash

for i in *.cssr;
do
  mkdir ${i%.cssr}
  mv $i ${i%.cssr}/.
  cp *.py ${i%.cssr}/.
  cp interp-generator-submit.sh ${i%.cssr}/.
  cd ${i%.cssr}
  qsub -N ${i%.cssr} interp-generator-submit.sh
  sleep 1m
  cd ../
done
