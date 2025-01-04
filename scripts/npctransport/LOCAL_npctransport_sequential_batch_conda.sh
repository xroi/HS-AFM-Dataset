#!/bin/bash
#SBATCH --mem=2g
#SBATCH --time=14-0
#SBATCH --array=1-25
#SBATCH --killable
#SBATCH --requeue
#SBATCH -c2
# Make sure to give array size with --array=1-x/2 (for example: --array=1-50 for 100 simulations)

if [ $# != 8 ]; then
  echo "Syntax: $0 <from-start(1/0)> <start> <step> <output_folder_path> <config_path> <output_statistics_interval> <stop_at(excluding)> <divergence(1/0)>"
  exit 1
fi

SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/

mkdir -p ${4}

declare -a IDs=()
IDs+=($((${SLURM_ARRAY_TASK_ID} * 2 - 1)))
IDs+=($((${SLURM_ARRAY_TASK_ID} * 2)))
echo IDs: ${IDs[@]}

echo "Running jobs"
declare -a PIDs=()
for ID in ${IDs[@]}; do
  mkdir -p ${4}/${ID}
  # If it's a divergence, copy the config file to the folder with correct timestamp naming.
  if [ "${8}" -eq "1" ]; then
        cp ${5} ${4}/${ID}/${2}.pb
  fi
  ${SCRIPTS_FOLDER}/npctransport/LOCAL_npctransport_sequential_conda.sh ${1} ${2} ${3} ${4}/${ID} ${5} ${6} ${7} &
  PID=$!
  PIDs+=($PID)
  echo $PID submitted, workid $ID
done

echo "Waiting for all jobs to finish"
for PID in "${PIDs[@]}"; do
  wait "$PID"
  echo "$PID finished"
done