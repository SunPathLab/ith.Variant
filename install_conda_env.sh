#!/bin/bash
set -e  

CODE_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing ith.Variant conda env..."

conda create -n ith.variant --file ${CODE_DIR}/requirements.txt \
  --override-channels -c bioconda -c conda-forge -c defaults -y

source activate ith.variant

echo "$(date): Installing modified version of TitanCNA..."

R -e 'install.packages("pkgs/TitanCNA_1.26.0.tar.gz")'

echo "$(date): Downloading muTect v1.1.5 jar ./pkgs folder..."

## install mutect v1.1.5
cd pkgs && wget --quiet https://github.com/broadinstitute/mutect/releases/download/1.1.5/muTect-1.1.5-bin.zip \
 && unzip muTect-1.1.5-bin.zip \
 && rm muTect-1.1.5-bin.zip && cd ..

echo "$(date): Installed ith.Variant Successfully!!"
