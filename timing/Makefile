NVCC=nvcc
LDLIBS=-lcusparse -lcusolver
DATASET_DIR=./baremetal_tests/incremental_sphere2500_steps-2-200_period-25

# List all header files in DATASET_DIR
HEADERS=$(wildcard $(DATASET_DIR)/*.h)

test_multithreaded_incremental_dataset_no_values_gpu: test_multithreaded_incremental_dataset_no_values_gpu.cu
	$(NVCC) -o $@ $< $(HEADERS) $(LDFLAGS)
