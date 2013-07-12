#!/bin/bash

cd build
make clean
make -j 2
if [ $? != 0 ]; then
  echo "ERROROROR"
  cd - &> /dev/null
  exit
fi
make install
cd - &> /dev/null

#./bin/vdw -a 1 -d 0.5 test.cube
#./bin/vdw -a 3 -d 0.5 bq1-homo0-100.cube
./bin/vdw -c -a 3 -d 1.0 bq1-homo0-100.cube
