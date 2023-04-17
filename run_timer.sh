VAR1="reordered_parking-garage"
VAR2="reordered_sphere_bignoise_vertex3"
# VAR3="reordered_grid3D"
# VAR4="w10000.graph"

read -p "Num threads (int): " NUM_THREADS
# NUM_THREADS=16

echo "Timing $VAR1"
./timing/testGtsamIncremental3D -f $VAR1 -k 1 > tbb-$NUM_THREADS-incremental-$VAR1.out

echo "Timing $VAR2"
./timing/testGtsamIncremental3D -f $VAR2 -k 1 > tbb-$NUM_THREADS-incremental-$VAR2.out

# echo "Timing $VAR3"
# ./timing/testGtsamIncremental3D -f $VAR3 -k 1 > no_tbb_incremental-$VAR3.out

# echo "Timing $VAR4"
# ./timing/testGtsamIncremental3D -f $VAR4 -k 1 > no_tbb_incremental-$VAR4.out
