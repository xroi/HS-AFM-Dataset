#!/bin/bash

if [ $# != 1 ]; then
  echo "Syntax: $0 <json_path>"
  exit 1
fi

SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/

length=$(jq length ${1})

for ((i=0; i<length; i++))
do
  amount=$(jq '.[${i}].amount' ${1})
  input_path=$(jq '.[${i}].input_path' ${1})
  output_path=$(jq '.[${i}].output_path' ${1})
  config_path=$(jq '.[${i}].config_path' ${1})

  sbatch --array=1-$((${amount}/2)) ${SCRIPTS_FOLDER}/afm/LOCAL_afm_batch.sh ${input_path} ${output_path} ${config_path}
done