#!/bin/bash
set -e  

CODE_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "Installing ith.Variant conda env...  [$(date)]"

conda create -n ith.variant --file ${CODE_DIR}/requirements.txt \
  --override-channels -c bioconda -c conda-forge -c defaults -y


echo "Activating conda env ith.variant...  [$(date)]"

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ith.variant

echo "Installing modified version of TitanCNA...  [$(date)]"

R -e 'install.packages("pkgs/TitanCNA_1.26.0.tar.gz")'

## install mutect v1.1.5
#echo "$(date): Downloading muTect v1.1.5 jar ${CODE_DIR}/pkgs folder..."
# cd pkgs && wget --quiet https://github.com/broadinstitute/mutect/releases/download/1.1.5/muTect-1.1.5-bin.zip \
#  && unzip muTect-1.1.5-bin.zip \
#  && rm muTect-1.1.5-bin.zip && cd ..

echo "ith.Variant installed successfully!! Please run: conda activate ith.variant"
