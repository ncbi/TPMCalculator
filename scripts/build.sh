#!/bin/bash

echo "Installing BAMtools"
git clone https://github.com/pezmaster31/bamtools
cd bamtools
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../ ..
make
make install
cd ..
export BAMTOOLS_DIR=`pwd`
export CPPFLAGS="-I $BAMTOOLS_DIR/include/bamtools/"
export LDFLAGS="-L $BAMTOOLS_DIR/lib -Wl,-rpath,$BAMTOOLS_DIR/lib"

echo "Using FLAGS:"
echo $BAMTOOLS_DIR
echo $CPPFLAGS
echo $LDFLAGS

echo "Compiling TPMCalculator"
pwd
ls
echo
echo
make


