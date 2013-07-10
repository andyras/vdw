#!/bin/bash

cd build
make -j 2 && make install
cd - &> /dev/null

./bin/vdw -a 0 -d 0.1 test.cube
