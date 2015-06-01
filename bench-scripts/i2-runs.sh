#!/bin/bash

CMD="/home/fernie/build/CMakeLists.txt-Desktop_GCC-Default/ih-generate-models blastOutput2.xml example2.fa"

for j in {1..10}; do
    { time $CMD; } 2>> i2-results${t}.txt;
done;

