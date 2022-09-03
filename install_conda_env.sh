#!/bin/bash
set -e  

CODE_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing ith.Variant conda env..."

conda create -n ith.variant --file ${CODE_DIR}/requirements.txt \
  --override-channels -c bioconda -c conda-forge -c defaults -c compbiocore -y

echo "$(date): Installed ith.Variant Successfully!!"
