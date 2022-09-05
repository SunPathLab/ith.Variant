#!/bin/bash
set -e  

CODE_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing ith.Variant conda env..."

conda create -n ith.variant --file ${CODE_DIR}/requirements.txt \
  --override-channels -c bioconda -c conda-forge -c defaults -c compbiocore -y

echo "$(date): Installed ith.Variant Successfully!!"

cd pkgs
zlib_version=1.2.12
wget http://www.zlib.net/zlib-${zlib_version}.tar.gz
tar -xvzf zlib-${zlib_version}.tar.gz

cd zlib-${zlib_version}

./configure --prefix=./zlib

make install

#wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz


