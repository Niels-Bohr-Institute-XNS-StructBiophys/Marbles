#!bin/bash

COLS=0
make
./runner > out_check

cd ../bead_modeling_flat/
gcc mainmc.c -o mainMC -lgsl -lgslcblas -L/usr/local/lib/ -I/usr/local/include/
./mainMC p450.input > out_check

cd ../c++/
echo ""
echo "RESULTS:"
python compare_outputs.py $COLS
