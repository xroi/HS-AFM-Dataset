#!/bin/csh -f

set IMP_FOLDER=/cs/labs/ravehb/ravehb/imp/fast_conda2/
set IMP=$IMP_FOLDER/setup_environment.sh
source /cs/labs/ravehb/ravehb/External/venv_imp2023_v2/bin/activate.csh

$IMP python3 $1 $2
