#!/bin/bash

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# for (( i=0; i<=40; i+=5 )); do
for (( i=0; i<=40; i+=5 )); do
    start_step=$i
    end_step=$(($start_step+4))
    if (( $start_step == 40 )); then
        end_step=44
    # if (( $start_step == 70 )); then
    #     end_step=73
    fi
    echo $(seq $start_step $end_step)
    python3 $SCRIPT_DIR/trim_header_nodes.py --infile $SCRIPT_DIR/../baremetal_tests/incremental_kitti_00_2LC_steps-863-863_period-1/step-863.h \
        --outfile $SCRIPT_DIR/../baremetal_tests/incremental_kitti_00_2LC_steps-863-863_period-1/step-863_node-${start_step}-${end_step}.h \
        --indices $(seq $start_step $end_step)
done

