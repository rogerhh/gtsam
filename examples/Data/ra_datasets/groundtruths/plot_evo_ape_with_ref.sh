#!/bin/bash

set -e

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

TRAJ_DIR=$1
REF_DIR=$2

traj_base=$(basename $TRAJ_DIR)
ref_base=$(basename $REF_DIR)

out_file=${traj_base}-${ref_base}_ape.out

traj_files=$(ls $TRAJ_DIR/step-*_traj.txt)

num_traj=$(ls -1 $TRAJ_DIR/step-*_traj.txt | wc -l)

step=0
endstep=3500

echo "Writing result to $out_file"
echo "" > $out_file

while [ $step -lt $endstep ]
do
  traj_file=$TRAJ_DIR/step-${step}_traj.txt
  ref_file=$REF_DIR/step-${step}_traj.txt

  echo $files_found 
  if [ -f ${traj_file} ]
  then

    ((files_found+=1))

    if [ -f ${ref_file} ]
    then
        echo $traj_file 
        echo $ref_file

        linecount=$(wc -l < $traj_file)
        tmp_ref=$1_$linecount
        head -n $linecount $ref_file > $tmp_ref
        evo_ape kitti $traj_file $tmp_ref >> $out_file

        rm $tmp_ref
    fi
  fi
  ((step+=1))
done

echo "Wrote result to $out_file"

