#!/bin/csh -f
#SBATCH --mem=4g
#SBATCH --time=7-0
#SBATCH --mail-type=ALL

if ($#argv != 3) then
    echo "Syntax: $0 <input_folder_path> <output_prefix> <args_path>"
    exit 0
endif

# Make sure you don't supply the --input-path and ----output-path-prefix arguments in the config file.

set AFM_SIM_PATH=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Simulation/

source $AFM_SIM_PATH/venv_new/bin/activate.csh

set INPUT_PATH=$1/
set OUTPUT_PATH=$2
set ORIGINAL_ARGS_PATH=$3
set ARGS_PATH=$OUTPUT_PATH._ARGS.txt
echo input path is $INPUT_PATH
echo output path is $OUTPUT_PATH

cp $ORIGINAL_ARGS_PATH $ARGS_PATH
echo "" >> $ARGS_PATH
echo "--input-path $INPUT_PATH" >> $ARGS_PATH
echo "--output-path-prefix $OUTPUT_PATH" >> $ARGS_PATH

python3 $AFM_SIM_PATH/src/main.py @$ARGS_PATH

