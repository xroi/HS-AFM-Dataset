#!/bin/csh -f
#SBATCH --mem=2g
#SBATCH --time=7-0
#SBATCH --killable
#SBATCH --requeue
# Make sure to give array size with --array=1-x/2 (for example: --array=1-50 for 100 simulations)

if ($#argv != 3) then
    echo "Syntax: $0 <input_folder_path> <output_folder_path> <config_path>"
    exit 0
endif

set INPUT_PATH=$1/
set OUTPUT_PATH=$2/
set CONFIG_PATH=$3
set SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/
# Note: Input path should be a folder, that has folder named 0-array_max
# (generated with LOCAL_npctransport_sequential_batch.sh)

mkdir $OUTPUT_PATH

mkdir -p $OUTPUT_PATH
$SCRIPTS_FOLDER/afm/LOCAL_afm.sh $INPUT_PATH/$SLURM_ARRAY_TASK_ID $OUTPUT_PATH/$SLURM_ARRAY_TASK_ID $CONFIG_PATH

