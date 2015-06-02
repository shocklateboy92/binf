#!/bin/bash

CMD="/home/fernie/build/CMakeLists.txt-Desktop_GCC-Default/ih-generate-models blastOutput2.xml example2.fa"

for t in {1..12}; do
    for j in {1..20}; do
        { time $CMD $t; } 2>> i2-results${t}.txt;
    done;
done;
