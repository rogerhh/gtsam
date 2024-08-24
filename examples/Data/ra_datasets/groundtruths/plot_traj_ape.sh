#!/bin/bash

set -e

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIO_CAB1-incremental_CAB1_ape.out VIOLC_CAB1-incremental_CAB1_ape.out RA_CAB1_num_threads-1-incremental_CAB1_ape.out RA_CAB1_num_threads-2-incremental_CAB1_ape.out RA_CAB1_num_threads-4-incremental_CAB1_ape.out " --labels="VIO VIO+LC RA-1thread RA-2thread RA-4thread" --outfile CAB464_traj_ape.png

# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIO_CAB1_vio_lag-20-LC_CAB1_ape.out \
#   VIOLC_CAB1_vio_lag-20_lc_period-8-LC_CAB1_ape.out RA_CAB1_num_threads-1-LC_CAB1_ape.out RA_CAB1_num_threads-2-LC_CAB1_ape.out RA_CAB1_num_threads-4-LC_CAB1_ape.out incremental_CAB1-LC_CAB1_ape.out " --labels="VIO VIO+LC RA-1thread RA-2thread RA-4thread Incremental" --outfile CAB464-LC_traj_ape.png
# 
# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIO_sphere2500_vio_lag-20-LC_sphere2500_ape.out VIOLC_sphere2500_vio_lag-20_lc_period-25-LC_sphere2500_ape.out RA_sphere2500_num_threads-1-LC_sphere2500_ape.out RA_sphere2500_num_threads-2-LC_sphere2500_ape.out RA_sphere2500_num_threads-4-LC_sphere2500_ape.out incremental_sphere2500-LC_sphere2500_ape.out" --labels="VIO VIO+LC RA-1thread RA-2thread RA-4thread Incremental" --outfile sphere2500-LC_traj_ape.png
# 
# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIO_CAB7000_vio_lag-20-LC_CAB7000_ape.out VIOLC_CAB7000_vio_lag-20_lc_period-65-LC_CAB7000_ape.out RA_CAB7000_num_threads-1-LC_CAB7000_ape.out RA_CAB7000_num_threads-2-LC_CAB7000_ape.out RA_CAB7000_num_threads-4-LC_CAB7000_ape.out incremental_CAB7000-LC_CAB7000_ape.out" --labels="VIO VIO+LC RA-1thread RA-2thread RA-4thread Incremental" --outfile CAB7000-LC_traj_ape.png
# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIOLC_CAB7000_vio_lag-20_lc_period-65-LC_CAB7000_ape.out RA_CAB7000_num_threads-4-LC_CAB7000_ape.out" --labels="VIO+LC RA-4thread" --outfile CAB7000-LC_traj_ape.png
# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIOLC_CAB7000_vio_lag-20_lc_period-65-incremental_CAB7000_ape1.out \
#   RA_CAB7000_num_threads-1-incremental_CAB7000_ape1.out \
#   RA_CAB7000_num_threads-2-incremental_CAB7000_ape1.out \
#   RA_CAB7000_num_threads-4-incremental_CAB7000_ape1.out" --labels="VIO+LC RA-1thread RA-2thread RA-4thread" --outfile CAB7000-traj_ape1.png

# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "VIO_CSAIL-LC_CSAIL_ape.out VIOLC_CSAIL_vio_lag-20_lc_period-65-LC_CSAIL_ape.out \
#   RA_CSAIL_num_threads-1-LC_CSAIL_ape.out \
#   RA_CSAIL_num_threads-2-LC_CSAIL_ape.out \
#   RA_CSAIL_num_threads-4-LC_CSAIL_ape.out" --labels="VIO VIO+LC RA-1thread RA-2thread RA-4thread" --outfile CSAIL-traj_ape1.png

  # VIO_CAB7000-smallnoise_vio_lag-20-LC_CAB7000-smallnoise_ape.out \
    # VIO_sphere2_vio_lag-20-LC_sphere2_ape.out \
  # VIO_CAB1-smallnoise_vio_lag-20-LC_CAB1-smallnoise_ape.out \
  # VIO_CAB7000-smallnoise_vio_lag-20-LC_CAB7000-smallnoise_ape.out \

python3 $SCRIPT_DIR/plot_traj_ape0.py --infiles "
    VIO_sphere2_vio_lag-20-LC_sphere2_ape.out \
    VIOLC_sphere2_vio_lag-20_lc_period-65-LC_sphere2_ape.out \
    RA_sphere2_num_threads-1-LC_sphere2_ape.out \
    RA_sphere2_num_threads-2-LC_sphere2_ape.out \
    RA_sphere2_num_threads-4-LC_sphere2_ape.out \
    incremental_sphere2-LC_sphere2_ape.out" \
    --labels="VIO VIO+LC RA1S RA2S RA4S In" --outfile sphere2-LC_traj_ape.png

# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "
#   VIOLC_CAB1-smallnoise_vio_lag-20_lc_period-8-LC_CAB1-smallnoise_ape.out \
#   RA_CAB1-smallnoise_num_threads-1-LC_CAB1-smallnoise_ape.out \
#   RA_CAB1-smallnoise_num_threads-2-LC_CAB1-smallnoise_ape.out \
#   RA_CAB1-smallnoise_num_threads-4-LC_CAB1-smallnoise_ape.out \
#   incremental_CAB1-smallnoise-LC_CAB1-smallnoise_ape.out" \
#     --labels="VIO+LC RA-1thread RA-2thread RA-4thread Incremental" --outfile CAB1-smallnoise-traj_ape.png

# python3 $SCRIPT_DIR/plot_traj_ape0.py --infiles "
#   VIO_CAB7000-smallnoise_vio_lag-20-LC_CAB7000-smallnoise_ape.out \
#   VIOLC_CAB7000-smallnoise_vio_lag-20_lc_period-65-LC_CAB7000-smallnoise_ape.out \
#   RA_CAB7000-smallnoise_num_threads-1-LC_CAB7000-smallnoise_ape.out \
#   RA_CAB7000-smallnoise_num_threads-2-LC_CAB7000-smallnoise_ape.out \
#   RA_CAB7000-smallnoise_num_threads-4-LC_CAB7000-smallnoise_ape.out \
#   incremental_CAB7000-smallnoise-LC_CAB7000-smallnoise_ape.out" \
#     --labels="VIO VIO+LC RA1S RA2S RA4S In" --outfile CAB7000-smallnoise-traj_ape.png
# 
# python3 $SCRIPT_DIR/plot_traj_ape.py --infiles "
#   LRU_CAB7000-smallnoise_num_threads-2_memsize-2M-LC_CAB7000-smallnoise_ape.out \
#   LRU_CAB7000-smallnoise_num_threads-2_memsize-4M-LC_CAB7000-smallnoise_ape.out \
#   VIO_CAB7000-smallnoise_vio_lag-716-LC_CAB7000-smallnoise_ape.out \
#   VIO_CAB7000-smallnoise_vio_lag-1024-LC_CAB7000-smallnoise_ape.out" \
#     --labels="LRU2M LRU4M VIO2M VIO4M" --outfile CAB7000-smallnoise-LRU-traj_ape.png

# python3 $SCRIPT_DIR/plot_traj_ape1.py --infiles "
#   VIO_M3500a_vio_lag-20-LC_M3500a_ape.out \
#   VIOLC_M3500a_vio_lag-20_lc_period-65-LC_M3500a_ape.out \
#   RA_M3500a_num_threads-1-LC_M3500a_ape.out \
#   RA_M3500a_num_threads-2-LC_M3500a_ape.out \
#   RA_M3500a_num_threads-4-LC_M3500a_ape.out \
#   incremental_M3500a-LC_M3500a_ape.out" \
#     --labels="VIO VIO+LC RA1S RA2S RA4S In" --outfile M3500a-traj_ape.png

python3 $SCRIPT_DIR/plot_traj_ape1.py --infiles "
  VIO_M3500a_vio_lag-20-LC_M3500a_ape.out \
  VIOLC_M3500a_vio_lag-20_lc_period-65-LC_M3500a_ape.out \
  RA_M3500a_num_threads-1-LC_M3500a_ape.out \
  RA_M3500a_num_threads-2-LC_M3500a_ape.out \
  RA_M3500a_num_threads-4-LC_M3500a_ape.out \
  incremental_M3500a-LC_M3500a_ape.out" \
    --labels="VIO VIO+LC RA1S RA2S RA4S In" --outfile M3500a-traj_ape.png

