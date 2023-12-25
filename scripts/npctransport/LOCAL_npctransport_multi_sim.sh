#!/bin/bash

if [ $# != 1 ]; then
  echo "Syntax: $0 <json_path>"
  exit 1
fi

SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/

length=$(jq length ${1})

for ((i=0; i<length; i++))
do
  amount=$(jq ".[${i}].amount" ${1})
  config_path=$(jq ".[${i}].config_path" ${1})
  output_path=$(jq ".[${i}].output_path" ${1})
  step=$(jq ".[${i}].step" ${1})
  stop=$(jq ".[${i}].stop" ${1})
  output_statistics_interval=$(jq ".[${i}].output_statistics_interval" ${1})

  sbatch --array=1-$((${amount}/2)) ${SCRIPTS_FOLDER}/npctransport/LOCAL_npctransport_sequential_batch.sh ${output_path} ${config_path} ${step} ${output_statistics_interval} ${stop}
done