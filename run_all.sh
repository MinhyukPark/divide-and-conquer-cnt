#!/bin/bash

/usr/bin/time -v python ./src/main.py --input-cnp-dir ./input/ --output-prefix ./output/k4_n15/ 2> errors/k4_n15/main.err 1> output/k4_n15/main.out

