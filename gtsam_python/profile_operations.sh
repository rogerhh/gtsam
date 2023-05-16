#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
DATA_DIR="${SCRIPT_DIR}/data"

datasets="sphere w10000"
opt="inc bat"

for d in $datasets
do
    for o in $opt
    do 
        python3 profile_operations.py --input_file1 "$DATA_DIR"/"$d"_"$o"_hotspots.csv --input_file2 "$DATA_DIR"/"$d"_"$o"_hw-events.csv --output_file "$DATA_DIR"/"$d"_"$o"_linalg-profile.csv
    done
done
