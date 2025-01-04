#!/bin/csh -f
#SBATCH --mem=2g
#SBATCH --time=7-0
#SBATCH --mail-type=ALL

# ===========
# Parameters:
# ===========
# $1 (from-start): if 1, initializes a new simulation. In this case, start should be equal to step. 
# $2 (start): used if from-start is 0. start.pb should be the the final pb file in output_folder_path.
# $3 (step): how often to restart simulations.
# $4 (output_folder_path): path of output folder.
# $5 (config_path): path of config (.pb) file ro be used for the simulations.
# $6 (output_statistics_interval): should be the same as in the pb file. For cleaning hdf5 names.
# $7 (stop_at): stops the simulation at this time (excluding).
#
# Important:
# - If from-start is 0, config_path isn't used. 
# - If from-start is 0, and start is -1, determines the start time dynamically by looking at the files in output_folder_path.
# - If from-start is 0, and start is i>0, the file output_folder_path/i.pb is assumed to exist.
# ===========

if ($#argv != 7) then
    echo "Syntax: $0 <from-start(1/0)> <start> <step> <output_folder_path> <config_path> <output_statistics_interval> <stop_at>"
    exit 0
endif

# Setup Environment
source /cs/labs/ravehb/roi.eliasian/miniconda3/etc/profile.d/conda.csh
conda activate imp_conda

set OUTPUT_PATH=$4/
set CONFIG_PATH=$5
set seed=`od -An -N4 -td4 /dev/random`
set SCRIPTS_FOLDER=/cs/labs/ravehb/roi.eliasian/NpcTransportExperiment/HS-AFM-Dataset/scripts/
set IMP_FOLDER=/cs/labs/ravehb/roi.eliasian/miniconda3/envs/imp_conda/

set start=$2

mkdir -p $OUTPUT_PATH
echo output path is $OUTPUT_PATH

if (`echo "$1==1" | bc`) then
    echo Initilising new simulation...
    set i=$3
    $IMP_FOLDER/bin/fg_simulation --configuration $CONFIG_PATH --output $OUTPUT_PATH$i.pb --short_init_factor 0.5 --short_sim_factor 1.00 --conformations $OUTPUT_PATH$i.movie.rmf --final_conformations $OUTPUT_PATH$i.pb.final.rmf --random_seed $seed
    ${SCRIPTS_FOLDER}/npctransport/LOCAL_change_hdf5_names.sh $OUTPUT_PATH $i $3 $6
endif

if (`echo "$2==-1" | bc`) then
    # determines the start time by looking at the files in output_folder_path
    # get filename of second largest pb (if exists, otherwise get largest) 
    set fname=`tree -fQFi --sort=size $4 | grep .pb\" | head -2 | tail -1`
    # extract time from file name
    set start=`echo $fname | sed -E 's|.*/([0-9]+)\.pb"?|\1|'`
    set 
    # delete files 1 over
    set one_over=`echo "$start + $3" | bc`
    rm $4/$one_over.pb
    rm $4/$one_over.pb.hdf5
    rm $4/$one_over.pb.final.rmf
    rm $4/$one_over.movie.rmf
endif

set i=(`echo "$start + $3" | bc`)
set j=$start
while (`echo "$i!=$7" | bc`)
    $IMP_FOLDER/bin/fg_simulation --output $OUTPUT_PATH$i.pb --conformations $OUTPUT_PATH$i.movie.rmf --final_conformations $OUTPUT_PATH$i.pb.final.rmf --restart $OUTPUT_PATH$j.pb
    ${SCRIPTS_FOLDER}/npctransport/LOCAL_keep_biggest_pb.sh $OUTPUT_PATH
    ${SCRIPTS_FOLDER}/npctransport/LOCAL_change_hdf5_names.sh $OUTPUT_PATH $i $3 $6
    echo cur: $i using:$j
    @ i+=$3
    @ j+=$3
end
