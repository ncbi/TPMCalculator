#!/bin/bash

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
export LDFLAGS="-L $BAMTOOLS_DIR/lib64 -Wl,-rpath,$BAMTOOLS_DIR/lib64"
make


