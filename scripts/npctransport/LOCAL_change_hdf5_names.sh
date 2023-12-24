#!/bin/bash

if [ $# != 4 ]; then
  echo "Syntax: $0 <dir_path> <time> <step> <output_statistics_interval>"
  exit 1
fi

time=${2}
step=${3}
interval=${4}
total=$((step/interval))

for (( i=1 ; i<=${total} ; i++ ));
do
  cur=$(((time - step) + (i * interval)))
  mv ${1}/${time}.pb.${i}.hdf5 ${1}/${cur}.hdf5
done