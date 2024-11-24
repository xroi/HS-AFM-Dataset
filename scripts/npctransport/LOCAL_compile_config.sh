#!/bin/csh -f
source ~/.zshrc
conda activate imp_conda

python3 $1 $2
