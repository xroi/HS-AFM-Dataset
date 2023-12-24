#!/bin/bash
#SBATCH --mem=2g
#SBATCH --time=7-0
#SBATCH --array=1-50
#SBATCH --killable
#SBATCH --requeue
#SBATCH -c2
# Make sure to give array size with --array=1-x/2 (for example: --array=1-50 for 100 simulations)

SECOND_JOB_OFFSET=50

if [ $# != 5 ]; then
  echo "Syntax: $0 <output_folder_path> <config_path> <step> <output_statistics_interval> <stop>"
  exit 1
fi

OUTPUT_PATH=${1}/
CONFIG_PATH=${2}
STEP=${3}
SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/

mkdir -p $OUTPUT_PATH

declare -a IDs=()
IDs+=(${SLURM_ARRAY_TASK_ID})
IDs+=($((SLURM_ARRAY_TASK_ID + SECOND_JOB_OFFSET)))
echo IDs: ${IDs[@]}

echo "Running jobs"
declare -a PIDs=()
for ID in ${IDs[@]}; do
  ${SCRIPTS_FOLDER}/npctransport/LOCAL_npctransport_sequential.sh 1 $STEP $STEP $OUTPUT_PATH/${ID} $CONFIG_PATH $4 $5 &
  PID=$!
  PIDs+=($PID)
  echo $PID submitted, workid $ID
done

echo "Waiting for all jobs to finish"
for PID in "${PIDs[@]}"; do
  wait "$PID"
  echo "$PID finished"
done