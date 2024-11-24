#!/usr/bin/env sh


LD_LIBRARY_PATH="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/lib:/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

PYTHONPATH="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/lib:/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/lib:$PYTHONPATH"
export PYTHONPATH

# Where to find data for the various modules
IMP_DATA="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/data"
export IMP_DATA

# Extra places to look for imp modules
IMP_EXAMPLE_DATA="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/doc/examples"
export IMP_EXAMPLE_DATA

# Location of binaries (for wine builds, which don't get PATH)
IMP_BIN_DIR="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/bin"
export IMP_BIN_DIR

PATH="/cs/labs/ravehb/ravehb/imp/fast_conda2/bin:/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/module_bin/algebra:/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/module_bin/atom:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/bayesianem:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/benchmark:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/bff:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/container:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/core:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/display:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/em:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/gsl:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/isd:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/kernel:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/kinematics:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/kmeans:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/misc:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/mmcif:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/modeller:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/mpi:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/multi_state:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/npc:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/npctransport:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/parallel:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/pmi:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/pmi1:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/rmf:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/rotamer:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/sampcon:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/saxs:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/saxs_merge:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/score_functor:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/scratch:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/statistics:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/symmetry:/cs/labs/ravehb/ravehb/imp/fast_conda2/module_bin/test:$PATH"
export PATH


IMP_TMP_DIR="/cs/labs/ravehb/ravehb/imp/Old/fast_conda2/tmp"
export IMP_TMP_DIR


mkdir -p "${IMP_TMP_DIR}"

exec ${precommand} "$@"