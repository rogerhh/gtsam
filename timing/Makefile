NVCC=nvcc
LDLIBS=-lcusparse -lcusolver
DATASET_DIR=./baremetal_tests/incremental_sphere2500_steps-2-200_period-25

test_multithreaded_incremental_dataset_no_values_gpu: test_multithreaded_incremental_dataset_no_values_gpu.cu
	$(NVCC) -o $@ $< -I$(DATASET_DIR) $(LDFLAGS)
