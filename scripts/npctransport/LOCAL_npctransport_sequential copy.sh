#!/bin/csh -f
#SBATCH --mem=2g
#SBATCH --time=7-0
#SBATCH --mail-type=ALL

if ($#argv != 7) then
    echo "Syntax: $0 <from-start(1/0)> <start> <step> <output_folder_path> <config_path> <output_statistics_interval> <stop_at(excluding)>"
    exit 0
endif

# Setup Environment
set CONDA_FOLDER=/cs/labs/ravehb/inbal/miniconda3/
source $CONDA_FOLDER/etc/profile.d/conda.csh
conda activate imp_conda

set OUTPUT_PATH=$4/
set CONFIG_PATH=$5
set seed=`od -An -N4 -td4 /dev/random`
# set SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/


mkdir -p $OUTPUT_PATH
echo output path is $OUTPUT_PATH

if (`echo "$1==1" | bc`) then
    echo Initilising new simulation...
    set i=$3
    $CONDA_FOLDER/envs/imp/bin/fg_simulation --configuration $CONFIG_PATH --output $OUTPUT_PATH$i.pb --short_init_factor 0.5 --short_sim_factor 1.00 --conformations $OUTPUT_PATH$i.movie.rmf --final_conformations $OUTPUT_PATH$i.pb.final.rmf --random_seed $seed
    # ${SCRIPTS_FOLDER}/npctransport/LOCAL_change_hdf5_names.sh $OUTPUT_PATH $i $3 $6
endif
