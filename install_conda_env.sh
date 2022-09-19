#!/bin/bash
set -e  

CODE_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing ith.Variant conda env..."

conda create -n ith.variant --file ${CODE_DIR}/requirements.txt \
  --override-channels -c bioconda -c conda-forge -c defaults -c compbiocore -y

source activate ith.variant

R -e 'install.packages("pkgs/TitanCNA_1.26.0.tar.gz")'

#conda install -c compbiocore mutect=1.1.6

echo "$(date): Installed ith.Variant Successfully!!"
