#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_DIR="$PROJECT_DIR/timing/python"
CONFIG_DIR="$PYTHON_DIR/configs"
BUILD_DIR="$PROJECT_DIR/build"

# python run_generate_dataset.py -y

logs=$(ls $BUILD_DIR/*.log)

cd $BUILD_DIR

# run=RA_sphere2_num_threads-1
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/sphere_smallnoise_pruned --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 1 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

# run=RA_sphere2_num_threads-2
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/sphere_smallnoise_pruned --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 2 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_sphere2_num_threads-4
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/sphere_smallnoise_pruned --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 4 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

# run=RA_CAB7000-smallnoise_num_threads-1
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB7000-2 --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 1 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_CAB7000-smallnoise_num_threads-2
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB7000-2 --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 2 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_CAB7000-smallnoise_num_threads-4
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB7000-2 --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 4 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

# run=RA_M3500a_num_threads-1
# ./timing/testGtsamIncremental-ra -f ra_datasets/M3500a --relin_thresh 100 --num_steps 3499 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 1 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_M3500a_num_threads-2
# ./timing/testGtsamIncremental-ra -f ra_datasets/M3500a --relin_thresh 100 --num_steps 3499 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 2 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_M3500a_num_threads-4
# ./timing/testGtsamIncremental-ra -f ra_datasets/M3500a --relin_thresh 100 --num_steps 3499 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 4 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

# run=RA_CAB1-smallnoise_num_threads-1
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB_query_val_phone-2 --relin_thresh 100 --num_steps 464 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 1 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_CAB1-smallnoise_num_threads-2
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB_query_val_phone-2 --relin_thresh 100 --num_steps 464 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 2 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
# 
# run=RA_CAB1-smallnoise_num_threads-4
# ./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB_query_val_phone-2 --relin_thresh 100 --num_steps 464 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --num_threads 4 --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

run=incremental_sphere2
./timing/testGtsamIncremental3D-ra -f ra_datasets/sphere_smallnoise_pruned --relin_thresh 100 --num_steps 2000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

run=incremental_CAB7000-smallnoise
./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB7000-2 --relin_thresh 100 --num_steps 3000 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

run=incremental_M3500a
./timing/testGtsamIncremental-ra -f ra_datasets/M3500a --relin_thresh 100 --num_steps 3499 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out

run=incremental_CAB1-smallnoise
./timing/testGtsamIncremental3D-ra -f ra_datasets/CAB_query_val_phone-2 --relin_thresh 100 --num_steps 464 --inject_delta_dir $run --skip_steps_file $run/skip_steps.txt --ra_latency_ms 33.3 --relin_keys_file $run/relin_keys.txt | tee $run.sym.out
