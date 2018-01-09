#!/bin/bash

# Command line arguments specify the sizes of the matrice under test.
# Possible values: 0008, 0016, 0064, 0256, 0512, 1024
#

# Specify where the input data are to be found, f.ex.
DATA_DIR=$(pwd)/data_in

while [ "$1" != "" ]; do
    SIZE=$1; shift 1
    echo "Testing ${DATA_DIR}/matrix_${SIZE}.mm"
    export OMP_NUM_THREADS=1
    echo "OMP_NUM_THREADS=1"
    step_5/step_5.x \
        -i ${DATA_DIR}/matrix_${SIZE}.mm -o step_5/res_${SIZE}_step_5_1.mm

    export OMP_NUM_THREADS=4
    echo "OMP_NUM_THREADS=4"
    step_5/step_5.x \
        -i ${DATA_DIR}/matrix_${SIZE}.mm -o step_5/res_${SIZE}_step_5_4.mm

    export OMP_NUM_THREADS=8
    echo "OMP_NUM_THREADS=8"
    step_5/step_5.x \
        -i ${DATA_DIR}/matrix_${SIZE}.mm -o step_5/res_${SIZE}_step_5_8.mm

    export OMP_NUM_THREADS=16
    echo "OMP_NUM_THREADS=16"
    step_5/step_5.x \
        -i ${DATA_DIR}/matrix_${SIZE}.mm -o step_5/res_${SIZE}_step_5_16.mm

    export OMP_NUM_THREADS=24
    echo "OMP_NUM_THREADS=24"
    step_5/step_5.x \
        -i ${DATA_DIR}/matrix_${SIZE}.mm -o step_5/res_${SIZE}_step_5_24.mm
    echo
done
