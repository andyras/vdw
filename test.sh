#!/bin/bash

cd build
make clean
make -j 2
if [ $? != 0 ]; then
  echo "ERROROROR"
  exit
fi
make install
cd - &> /dev/null

./bin/vdw -a 0 -d 0.1 test.cube
